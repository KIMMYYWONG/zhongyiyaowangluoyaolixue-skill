if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop(paste(
    "Usage: Rscript 06_build_network_tables.R",
    "<compound_target.csv|tsv> <intersection_targets.csv|tsv|txt>",
    "<formula_name> <disease_name> <herb_code_map.csv|tsv>",
    "<complete|filtered> <output_dir> [full_herb_compound_membership|NONE]"
  ))
}
options(stringsAsFactors = FALSE)
compound_target_file <- args[[1]]
intersection_file <- args[[2]]
formula_name <- trimws(args[[3]])
disease_name <- trimws(args[[4]])
herb_code_file <- args[[5]]
assignment_mode <- tolower(trimws(args[[6]]))
output_dir <- args[[7]]
membership_file <- if (length(args) >= 8) args[[8]] else "NONE"
if (!(assignment_mode %in% c("complete", "filtered"))) stop("assignment mode must be complete or filtered.")
if (!nzchar(formula_name) || !nzchar(disease_name)) stop("Formula and disease names must not be empty.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

split_input_spec <- function(spec) {
  hit <- regexpr("::[^:]+$", spec)
  if (hit[[1]] < 0) return(list(path = spec, sheet = NULL))
  list(path = substr(spec, 1, hit[[1]] - 1), sheet = substr(spec, hit[[1]] + 2, nchar(spec)))
}
read_table <- function(spec) {
  input <- split_input_spec(spec); path <- input$path
  if (!file.exists(path)) stop("Input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (ext %in% c("tsv", "txt")) {
    first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
    if (ext == "txt" && length(first) && !grepl("\t", first, fixed = TRUE)) return(data.frame(Gene_symbol = readLines(path, warn = FALSE, encoding = "UTF-8")))
    return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  }
  if (ext %in% c("xlsx", "xlsm")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Reading Excel input requires the openxlsx package.")
    sheets <- openxlsx::getSheetNames(path)
    if (is.null(input$sheet) && length(sheets) != 1L) {
      stop("Excel workbook has multiple sheets. Select one explicitly as file.xlsx::SheetName")
    }
    sheet <- if (is.null(input$sheet)) 1L else if (grepl("^[0-9]+$", input$sheet)) as.integer(input$sheet) else input$sheet
    return(openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE, detectDates = FALSE))
  }
  stop("Supported inputs: CSV, TSV, TXT, XLSX, and XLSM: ", path)
}
write_out <- function(x, name) write.csv(x, file.path(output_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
key <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
find_col <- function(df, candidates, required = TRUE) {
  hit <- match(key(candidates), key(names(df)), nomatch = 0); hit <- hit[hit > 0]
  if (!length(hit)) { if (required) stop("Missing column; expected one of: ", paste(candidates, collapse = ", ")); return(NA_character_) }
  names(df)[hit[[1]]]
}
clean_text <- function(x) {
  y <- trimws(gsub("[[:space:]]+", " ", as.character(x)))
  y[is.na(y)] <- ""
  y
}
value_or_blank <- function(df, col) if (is.na(col)) rep("", nrow(df)) else clean_text(df[[col]])

normalize_relations <- function(df, need_gene = TRUE) {
  herb_col <- find_col(df, c("Herb_name", "Herb", "Drug_name", "Medicine", "中药名"))
  gene_col <- find_col(df, c("Gene_symbol", "Gene Symbol", "Target_gene_symbol", "symbol", "target"), required = need_gene)
  name_col <- find_col(df, c("Compound_name", "Ingredient_name", "Ingredient", "Molecule_name", "成分名称", "Mol_name"))
  mol_col <- find_col(df, c("Mol_ID", "TCMSP_ID", "Ingredient_ID", "Compound_ID"), required = FALSE)
  pubchem_col <- find_col(df, c("PubChem_CID", "Pubchem Cid", "PubChem CID", "CID"), required = FALSE)
  inchikey_col <- find_col(df, c("InChIKey", "InChI Key"), required = FALSE)
  cas_col <- find_col(df, c("CAS", "CAS_number", "CAS No"), required = FALSE)
  herb <- clean_text(df[[herb_col]])
  gene <- if (need_gene) toupper(clean_text(df[[gene_col]])) else rep("", nrow(df))
  cname <- clean_text(df[[name_col]])
  mol <- value_or_blank(df, mol_col); pubchem <- value_or_blank(df, pubchem_col)
  inchikey <- toupper(value_or_blank(df, inchikey_col)); cas <- value_or_blank(df, cas_col)
  compound_key <- ifelse(inchikey != "", paste0("INCHIKEY:", inchikey),
    ifelse(pubchem != "", paste0("PUBCHEM:", pubchem),
      ifelse(cas != "", paste0("CAS:", cas), paste0("NAME:", key(cname)))))
  data.frame(Source_row = seq_len(nrow(df)), Herb_name = herb, Compound_key = compound_key,
    Compound_name = cname, Original_compound_ID = mol, PubChem_CID = pubchem,
    InChIKey = inchikey, CAS = cas, Gene_symbol = gene, stringsAsFactors = FALSE)
}

raw_rel <- read_table(compound_target_file)
rel_all <- normalize_relations(raw_rel, need_gene = TRUE)
rel_all <- rel_all[rel_all$Herb_name != "" & rel_all$Compound_name != "" & rel_all$Gene_symbol != "" & rel_all$Compound_key != "NAME:", ]
dedupe_columns <- setdiff(names(rel_all), "Source_row")
duplicate_flag <- duplicated(rel_all[, dedupe_columns, drop = FALSE])
exact_duplicates <- rel_all[duplicate_flag, , drop = FALSE]
rel <- rel_all[!duplicate_flag, , drop = FALSE]
write_out(exact_duplicates, "00_exact_duplicate_rows_removed.csv")
intersection <- read_table(intersection_file)
gene_col <- find_col(intersection, c("Gene_symbol", "Gene Symbol", "symbol", "gene", "preferredName", "name"))
intersection_genes <- sort(unique(toupper(clean_text(intersection[[gene_col]]))))
intersection_genes <- intersection_genes[intersection_genes != ""]
filtered <- rel[rel$Gene_symbol %in% intersection_genes, , drop = FALSE]
if (!nrow(filtered)) stop("No compound-target records match the frozen intersection genes.")

codes_raw <- read_table(herb_code_file)
code_herb_col <- find_col(codes_raw, c("Herb_name", "Herb", "Drug_name", "中药名"))
code_col <- find_col(codes_raw, c("Herb_code", "Code", "Abbreviation", "药物简称"))
codes <- unique(data.frame(Herb_name = clean_text(codes_raw[[code_herb_col]]), Herb_code = tolower(clean_text(codes_raw[[code_col]]))))
if (any(codes$Herb_name == "" | codes$Herb_code == "")) stop("Herb code map contains empty names or codes.")
if (anyDuplicated(codes$Herb_name) || anyDuplicated(codes$Herb_code)) stop("Herb names and herb codes must each be unique.")
missing_codes <- setdiff(unique(rel$Herb_name), codes$Herb_name)
if (length(missing_codes)) stop("Missing herb codes for: ", paste(missing_codes, collapse = ", "))

if (assignment_mode == "complete") {
  membership_raw <- if (toupper(membership_file) == "NONE") raw_rel else read_table(membership_file)
  membership <- normalize_relations(membership_raw, need_gene = FALSE)
  membership <- unique(membership[membership$Herb_name != "" & membership$Compound_key %in% filtered$Compound_key,
    c("Herb_name", "Compound_key", "Compound_name")])
} else {
  membership <- unique(filtered[, c("Herb_name", "Compound_key", "Compound_name")])
}
if (!nrow(membership)) stop("No herb-compound memberships remain under the selected assignment mode.")
missing_membership_codes <- setdiff(unique(membership$Herb_name), codes$Herb_name)
if (length(missing_membership_codes)) stop("Missing herb codes for membership herbs: ", paste(missing_membership_codes, collapse = ", "))

herb_lists <- aggregate(Herb_name ~ Compound_key, membership, function(x) paste(sort(unique(x)), collapse = ";"))
herb_counts <- aggregate(Herb_name ~ Compound_key, membership, function(x) length(unique(x)))
names(herb_counts)[[2]] <- "Herb_count"
ingredient_dict <- unique(filtered[, c("Compound_key", "Compound_name", "Original_compound_ID", "PubChem_CID", "InChIKey", "CAS")])
ingredient_dict <- ingredient_dict[order(ingredient_dict$Compound_key, ingredient_dict$Compound_name), ]
ingredient_dict <- ingredient_dict[!duplicated(ingredient_dict$Compound_key), ]
ingredient_dict <- merge(ingredient_dict, herb_lists, by = "Compound_key", all.x = TRUE)
ingredient_dict <- merge(ingredient_dict, herb_counts, by = "Compound_key", all.x = TRUE)
names(ingredient_dict)[names(ingredient_dict) == "Herb_name"] <- "herbs"
ingredient_dict$Shared <- ingredient_dict$Herb_count > 1
ingredient_dict$Component_node <- NA_character_
shared_idx <- which(ingredient_dict$Shared)
if (length(shared_idx)) ingredient_dict$Component_node[shared_idx] <- paste0("S", seq_along(shared_idx))
for (i in which(!ingredient_dict$Shared)) {
  herb <- ingredient_dict$herbs[[i]]
  hcode <- codes$Herb_code[match(herb, codes$Herb_name)]
  peers <- which(!ingredient_dict$Shared & ingredient_dict$herbs == herb)
  ingredient_dict$Component_node[[i]] <- paste0(hcode, match(i, peers))
}
ingredient_dict$Assignment_mode <- assignment_mode

filtered <- merge(filtered, ingredient_dict[, c("Compound_key", "Component_node", "herbs", "Shared")], by = "Compound_key", all.x = TRUE)
filtered <- filtered[order(filtered$Component_node, filtered$Compound_name, filtered$Gene_symbol, filtered$Herb_name), ]
filtered <- unique(filtered)
write_out(filtered, "01_number_ingredient_gene_audit.csv")
write_out(ingredient_dict, "02_ingredient_node_dictionary.csv")

formula_edges <- data.frame(source = formula_name, target = codes$Herb_name, interaction = "contains_drug")
membership_nodes <- merge(membership, ingredient_dict[, c("Compound_key", "Component_node")], by = "Compound_key")
herb_edges <- unique(data.frame(source = membership_nodes$Herb_name, target = membership_nodes$Component_node, interaction = "contains_ingredient"))
ingredient_edges <- unique(data.frame(source = filtered$Component_node, target = filtered$Gene_symbol, interaction = "targets"))
disease_edges <- data.frame(source = disease_name, target = intersection_genes, interaction = "associated_with")
edges <- unique(rbind(formula_edges, herb_edges, ingredient_edges, disease_edges))
edges <- edges[order(edges$interaction, edges$source, edges$target), ]

ingredient_nodes <- data.frame(term = ingredient_dict$Component_node, type = "Ingredient", label = ingredient_dict$Compound_name, herbs = ingredient_dict$herbs, entity_key = ingredient_dict$Compound_key)
nodes <- rbind(
  data.frame(term = formula_name, type = "Formula", label = formula_name, herbs = "", entity_key = paste0("FORMULA:", formula_name)),
  data.frame(term = codes$Herb_name, type = "Drug", label = codes$Herb_name, herbs = codes$Herb_name, entity_key = paste0("HERB:", codes$Herb_code)),
  ingredient_nodes,
  data.frame(term = intersection_genes, type = "Gene", label = intersection_genes, herbs = "", entity_key = paste0("GENE:", intersection_genes)),
  data.frame(term = disease_name, type = "Disease", label = disease_name, herbs = "", entity_key = paste0("DISEASE:", disease_name))
)

node_type_conflict <- aggregate(type ~ term, nodes, function(x) length(unique(x)))
conflicting_terms <- node_type_conflict$term[node_type_conflict$type > 1]
nodes <- nodes[!duplicated(nodes$term), ]
duplicate_edges <- anyDuplicated(edges[, c("source", "target", "interaction")]) > 0
blank_endpoints <- any(is.na(edges$source) | is.na(edges$target) | edges$source == "" | edges$target == "")
missing_endpoints <- setdiff(unique(c(edges$source, edges$target)), nodes$term)
shared_bad <- ingredient_dict$Component_node[(ingredient_dict$Shared & ingredient_dict$Herb_count < 2) | (!ingredient_dict$Shared & ingredient_dict$Herb_count != 1)]
expected_herb_edges <- unique(paste(membership_nodes$Herb_name, membership_nodes$Component_node, sep = "\r"))
actual_herb_edges <- unique(paste(herb_edges$source, herb_edges$target, sep = "\r"))
membership_mismatch <- length(setdiff(union(expected_herb_edges, actual_herb_edges), intersect(expected_herb_edges, actual_herb_edges)))

qc <- data.frame(Check = c("No blank edge endpoints", "No duplicate edges", "All endpoints in node table",
  "Unique node keys", "No conflicting node types", "Shared and unique ingredient herb counts valid", "Shared assignment matches selected mode"),
  Pass = c(!blank_endpoints, !duplicate_edges, !length(missing_endpoints), !anyDuplicated(nodes$term),
    !length(conflicting_terms), !length(shared_bad), membership_mismatch == 0),
  Detail = c("", "", paste(missing_endpoints, collapse = ";"), "", paste(conflicting_terms, collapse = ";"),
    paste(shared_bad, collapse = ";"), paste0("mode=", assignment_mode, "; mismatches=", membership_mismatch)))
counts <- data.frame(Stage = c("Raw compound-target rows", "Normalized valid rows", "Retained intersection rows",
  "Distinct herb-compound-gene", "Ingredient nodes", "Shared ingredient nodes", "Edges", "Nodes"),
  Count = c(nrow(raw_rel), nrow(rel), nrow(filtered), nrow(unique(filtered[, c("Herb_name", "Compound_key", "Gene_symbol")])),
    nrow(ingredient_dict), sum(ingredient_dict$Shared), nrow(edges), nrow(nodes)))
write_out(edges, "03_network_edges.csv")
write_out(nodes, "04_network_nodes.csv")
write_out(qc, "05_network_qc.csv")
write_out(counts, "06_combination_count_audit.csv")
write_out(membership_nodes[, c("Herb_name", "Component_node", "Compound_key")], "07_shared_assignment_audit.csv")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  workbook_tables <- list(Edges = edges, Nodes = nodes, Ingredient_dictionary = ingredient_dict, QC = qc, Counts = counts)
  for (nm in names(workbook_tables)) {
    openxlsx::addWorksheet(wb, nm); openxlsx::writeData(wb, nm, workbook_tables[[nm]])
  }
  openxlsx::saveWorkbook(wb, file.path(output_dir, "08_Cytoscape_network_tables.xlsx"), overwrite = TRUE)
}
if (!all(qc$Pass)) stop("Critical network QC failed; inspect 05_network_qc.csv before Cytoscape import.")
cat("assignment_mode=", assignment_mode, "\n", sep = "")
cat("ingredients=", nrow(ingredient_dict), "\n", sep = "")
cat("shared_ingredients=", sum(ingredient_dict$Shared), "\n", sep = "")
cat("edges=", nrow(edges), "\n", sep = "")
cat("nodes=", nrow(nodes), "\n", sep = "")
