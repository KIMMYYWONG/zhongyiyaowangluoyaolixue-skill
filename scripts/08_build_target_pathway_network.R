if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage: Rscript 08_build_target_pathway_network.R",
    "<ingredient_target.csv|xlsx|NONE> <target_set.csv|txt|xlsx>",
    "<KEGG_result.csv|xlsx> <intersection|core> <output_dir>",
    "[p_cutoff=0.05] [top_n=30]",
    "Use file.xlsx::SheetName for a multi-sheet workbook."
  ))
}
if (!requireNamespace("dplyr", quietly = TRUE)) stop("This script requires dplyr.")
if (!requireNamespace("openxlsx", quietly = TRUE)) stop("This script requires openxlsx for the Cytoscape workbook.")

ingredient_file <- args[[1]]
target_file <- args[[2]]
kegg_file <- args[[3]]
target_set <- tolower(trimws(args[[4]]))
output_dir <- args[[5]]
p_cutoff <- if (length(args) >= 6) as.numeric(args[[6]]) else 0.05
top_n <- if (length(args) >= 7) as.integer(args[[7]]) else 30L
if (!(target_set %in% c("intersection", "core"))) stop("target_set must be intersection or core.")
if (!is.finite(p_cutoff) || p_cutoff <= 0 || p_cutoff > 1) stop("p_cutoff must be in (0, 1].")
if (is.na(top_n) || top_n < 1L) stop("top_n must be a positive integer.")
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
    if (ext == "txt" && length(first) && !grepl("\t", first, fixed = TRUE)) {
      return(data.frame(Gene_symbol = readLines(path, warn = FALSE, encoding = "UTF-8")))
    }
    return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  }
  if (ext %in% c("xlsx", "xlsm")) {
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
  if (!length(hit)) {
    if (required) stop("Missing column; expected one of: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  names(df)[hit[[1]]]
}
clean <- function(x) {
  y <- trimws(gsub("[[:space:]]+", " ", as.character(x)))
  y[is.na(y)] <- ""
  y
}

targets_raw <- read_table(target_file)
target_col <- find_col(targets_raw, c("Gene_symbol", "Gene Symbol", "symbol", "gene", "preferredName", "name", "shared name"))
targets <- data.frame(Gene = unique(toupper(clean(targets_raw[[target_col]]))))
targets <- targets[targets$Gene != "", , drop = FALSE]
if (!nrow(targets)) stop("The selected target set is empty.")

kegg_raw <- read_table(kegg_file)
pathway_id_col <- find_col(kegg_raw, c("ID", "KEGGID", "Pathway_ID"))
pathway_name_col <- find_col(kegg_raw, c("Description", "Pathway", "Pathway_name"))
pathway_gene_col <- find_col(kegg_raw, c("Gene_symbols", "Gene_symbol", "geneID"))
p_col <- find_col(kegg_raw, c("pvalue", "P_value", "P.Value"))
padj_col <- find_col(kegg_raw, c("p.adjust", "p_adjust", "adjusted_pvalue", "FDR", "qvalue"), required = FALSE)

kegg <- data.frame(
  Pathway = clean(kegg_raw[[pathway_id_col]]),
  Pathway_name = clean(kegg_raw[[pathway_name_col]]),
  pvalue = suppressWarnings(as.numeric(kegg_raw[[p_col]])),
  p_adjust = if (is.na(padj_col)) NA_real_ else suppressWarnings(as.numeric(kegg_raw[[padj_col]])),
  Gene_list = clean(kegg_raw[[pathway_gene_col]]),
  stringsAsFactors = FALSE
)
kegg <- kegg[is.finite(kegg$pvalue) & kegg$pvalue < p_cutoff & kegg$Pathway != "" & kegg$Gene_list != "", , drop = FALSE]
kegg <- dplyr::slice_head(dplyr::arrange(kegg, .data$pvalue, .data$Pathway), n = top_n)
if (!nrow(kegg)) stop("No KEGG pathways pass the requested p-value cutoff.")

expanded <- lapply(seq_len(nrow(kegg)), function(i) {
  genes <- unique(toupper(trimws(unlist(strsplit(kegg$Gene_list[[i]], "[/;,]")))))
  genes <- genes[genes != ""]
  data.frame(Pathway = kegg$Pathway[[i]], Pathway_name = kegg$Pathway_name[[i]],
    pvalue = kegg$pvalue[[i]], p_adjust = kegg$p_adjust[[i]], Gene = genes,
    target_set = target_set, stringsAsFactors = FALSE)
})
pathway_gene_all <- dplyr::distinct(dplyr::bind_rows(expanded))
pathway_gene <- dplyr::inner_join(pathway_gene_all, targets, by = "Gene")
pathway_gene <- dplyr::arrange(dplyr::distinct(pathway_gene), .data$pvalue, .data$Pathway, .data$Gene)
if (!nrow(pathway_gene)) stop("No pathway genes match the explicitly selected target set.")

edges <- dplyr::distinct(dplyr::transmute(pathway_gene,
  source = .data$Pathway, target = .data$Gene, interaction = "contains_gene", target_set = .data$target_set))
pathway_nodes <- pathway_gene |>
  dplyr::arrange(.data$pvalue, .data$Pathway_name) |>
  dplyr::group_by(.data$Pathway, .data$target_set) |>
  dplyr::summarise(label = dplyr::first(.data$Pathway_name), pvalue = dplyr::first(.data$pvalue),
    p_adjust = dplyr::first(.data$p_adjust), .groups = "drop") |>
  dplyr::transmute(term = .data$Pathway, type = "Pathway", .data$label, .data$pvalue, .data$p_adjust, .data$target_set)
gene_nodes <- dplyr::distinct(dplyr::transmute(pathway_gene, term = .data$Gene, type = "Gene",
  label = .data$Gene, pvalue = NA_real_, p_adjust = NA_real_, target_set = .data$target_set))
nodes <- dplyr::bind_rows(pathway_nodes, gene_nodes)
nodes <- dplyr::arrange(dplyr::distinct(nodes), .data$type, .data$term)

triples <- data.frame()
if (toupper(ingredient_file) != "NONE") {
  ingredients_raw <- read_table(ingredient_file)
  ingredient_id_col <- find_col(ingredients_raw, c("Component_node", "Ingredient_node", "Mol_ID", "TCMSP_ID", "Ingredient_ID"))
  ingredient_name_col <- find_col(ingredients_raw, c("Compound_name", "Ingredient_name", "Ingredient", "Molecule_name", "成分名称", "Mol_name"), required = FALSE)
  ingredient_gene_col <- find_col(ingredients_raw, c("Gene_symbol", "Gene Symbol", "symbol", "target"))
  ingredient_gene <- data.frame(
    Ingredient = clean(ingredients_raw[[ingredient_id_col]]),
    Ingredient_name = if (is.na(ingredient_name_col)) clean(ingredients_raw[[ingredient_id_col]]) else clean(ingredients_raw[[ingredient_name_col]]),
    Gene = toupper(clean(ingredients_raw[[ingredient_gene_col]])), stringsAsFactors = FALSE
  )
  ingredient_gene <- dplyr::distinct(dplyr::filter(ingredient_gene, .data$Ingredient != "", .data$Gene != ""))
  triples <- dplyr::arrange(dplyr::distinct(dplyr::inner_join(
    ingredient_gene, pathway_gene, by = "Gene", relationship = "many-to-many"
  )),
    .data$Ingredient, .data$Gene, .data$Pathway)
}

blank_edges <- any(edges$source == "" | edges$target == "")
duplicate_edges <- anyDuplicated(edges[, c("source", "target", "interaction")]) > 0
missing_endpoints <- setdiff(unique(c(edges$source, edges$target)), nodes$term)
conflicting_types <- nodes |>
  dplyr::group_by(.data$term) |>
  dplyr::summarise(n_types = dplyr::n_distinct(.data$type), .groups = "drop") |>
  dplyr::filter(.data$n_types > 1L) |>
  dplyr::pull("term")
qc <- data.frame(
  Check = c("Target set explicitly labelled", "No blank edge endpoints", "No duplicate edges",
    "All endpoints in node table", "Unique node keys", "No conflicting node types",
    "KEGG genes restricted to selected target set"),
  Pass = c(target_set %in% c("intersection", "core"), !blank_edges, !duplicate_edges,
    !length(missing_endpoints), !anyDuplicated(nodes$term), !length(conflicting_types),
    all(pathway_gene$Gene %in% targets$Gene)),
  Detail = c(paste0("target_set=", target_set), "", "", paste(missing_endpoints, collapse = ";"),
    "", paste(conflicting_types, collapse = ";"), paste0("p<", p_cutoff, "; top_n=", top_n))
)
provenance <- data.frame(
  Item = c("target_set", "target_input", "KEGG_input", "ingredient_input", "p_cutoff", "top_n",
    "Excel_reader", "openxlsx_runtime_version", "join_engine"),
  Value = c(target_set, target_file, kegg_file, ingredient_file, p_cutoff, top_n,
    "openxlsx::read.xlsx when input is XLSX/XLSM", as.character(utils::packageVersion("openxlsx")), "dplyr")
)

write_out(edges, "01_target_pathway_edges.csv")
write_out(nodes, "02_target_pathway_nodes.csv")
write_out(pathway_gene, "03_pathway_gene_relations.csv")
if (nrow(triples)) write_out(triples, "04_ingredient_gene_pathway_triples.csv")
write_out(qc, "05_target_pathway_qc.csv")
write_out(provenance, "06_read_and_method_provenance.csv")

wb <- openxlsx::createWorkbook()
edge_legacy <- dplyr::transmute(edges, from_node = .data$source, to_node = .data$target)
node_legacy <- dplyr::transmute(nodes, node = .data$term, type = tolower(.data$type))
tables <- list(Edge_Table = edge_legacy, Node_Table = node_legacy,
  Canonical_Edges = edges, Canonical_Nodes = nodes, Pathway_Gene = pathway_gene,
  QC = qc, Provenance = provenance)
if (nrow(triples)) tables$Ingredient_Gene_Pathway <- triples
for (nm in names(tables)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, tables[[nm]])
}
openxlsx::saveWorkbook(wb, file.path(output_dir, "07_Cytoscape_target_pathway_tables.xlsx"), overwrite = TRUE)

if (!all(qc$Pass)) stop("Critical target-pathway QC failed; inspect 05_target_pathway_qc.csv.")
cat("target_set=", target_set, "\n", sep = "")
cat("pathways=", length(unique(pathway_gene$Pathway)), "\n", sep = "")
cat("genes=", length(unique(pathway_gene$Gene)), "\n", sep = "")
cat("edges=", nrow(edges), "\n", sep = "")
cat("triples=", nrow(triples), "\n", sep = "")
