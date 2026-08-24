if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 07_build_sankey.R <01_number_ingredient_gene_audit.csv> <KEGG_display.csv> <output_dir> [max_ingredients=30] [max_genes=30]")
relation_file <- args[[1]]; kegg_file <- args[[2]]; output_dir <- args[[3]]
max_ingredients <- if (length(args) >= 4) as.integer(args[[4]]) else 30L
max_genes <- if (length(args) >= 5) as.integer(args[[5]]) else 30L
if (max_ingredients < 1 || max_genes < 1) stop("Display limits must be positive integers.")
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
  if (ext %in% c("tsv", "txt")) return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  if (ext %in% c("xlsx", "xlsm")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Reading Excel input requires the openxlsx package.")
    sheets <- openxlsx::getSheetNames(path)
    if (is.null(input$sheet) && length(sheets) != 1L) {
      stop("Excel workbook has multiple sheets. Select one explicitly as file.xlsx::SheetName")
    }
    sheet <- if (is.null(input$sheet)) 1L else if (grepl("^[0-9]+$", input$sheet)) as.integer(input$sheet) else input$sheet
    return(openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE, detectDates = FALSE))
  }
  stop("Supported inputs: CSV, TSV, TXT, XLSX, and XLSM.")
}
write_out <- function(x, name) write.csv(x, file.path(output_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
key <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
find_col <- function(df, candidates) {
  hit <- match(key(candidates), key(names(df)), nomatch = 0); hit <- hit[hit > 0]
  if (!length(hit)) stop("Missing column; expected one of: ", paste(candidates, collapse = ", "))
  names(df)[hit[[1]]]
}
relations <- read_table(relation_file); kegg <- read_table(kegg_file)
component_col <- find_col(relations, c("Component_node", "Ingredient_node", "Mol_ID"))
component_name_col <- find_col(relations, c("Compound_name", "Ingredient_name", "Ingredient"))
gene_col <- find_col(relations, c("Gene_symbol", "Gene Symbol", "symbol"))
pathway_id_col <- find_col(kegg, c("ID", "KEGGID", "Pathway_ID"))
pathway_name_col <- find_col(kegg, c("Description", "Pathway", "Pathway_name"))
pathway_gene_col <- find_col(kegg, c("Gene_symbols", "Gene_symbol", "geneID"))
pathway_p_col <- find_col(kegg, c("pvalue", "P_value", "P.Value"))

ingredient_gene <- unique(data.frame(Ingredient = as.character(relations[[component_col]]),
  Ingredient_name = as.character(relations[[component_name_col]]), Gene = toupper(as.character(relations[[gene_col]]))))
pathway_gene_rows <- list()
for (i in seq_len(nrow(kegg))) {
  genes <- unlist(strsplit(as.character(kegg[[pathway_gene_col]][[i]]), "[/;,]"))
  genes <- toupper(trimws(genes)); genes <- genes[genes != ""]
  pathway_p <- suppressWarnings(as.numeric(kegg[[pathway_p_col]][[i]]))
  pathway_logp <- if (is.finite(pathway_p) && pathway_p > 0) -log10(pathway_p) else NA_real_
  if (length(genes)) pathway_gene_rows[[length(pathway_gene_rows) + 1L]] <- data.frame(
    Gene = genes, Pathway = as.character(kegg[[pathway_id_col]][[i]]),
    Pathway_name = as.character(kegg[[pathway_name_col]][[i]]),
    Pathway_pvalue = pathway_p, Pathway_minus_log10_p = pathway_logp,
    stringsAsFactors = FALSE)
}
if (!length(pathway_gene_rows)) stop("No pathway-gene relationships could be expanded from the KEGG display table.")
pathway_gene <- unique(do.call(rbind, pathway_gene_rows))
triples <- unique(merge(ingredient_gene, pathway_gene, by = "Gene"))
if (!nrow(triples)) stop("No ingredient-gene-pathway triples remain after matching.")
triples <- triples[order(triples$Ingredient, triples$Gene, triples$Pathway), ]
triples$weight <- 1
write_out(triples, "01_sankey_all_triples.csv")

ingredient_rank <- aggregate(Gene ~ Ingredient, triples, function(x) length(unique(x)))
names(ingredient_rank)[[2]] <- "Distinct_genes"; ingredient_rank <- ingredient_rank[order(-ingredient_rank$Distinct_genes, ingredient_rank$Ingredient), ]
keep_ingredients <- head(ingredient_rank$Ingredient, max_ingredients)
gene_rank <- aggregate(Pathway ~ Gene, triples[triples$Ingredient %in% keep_ingredients, ], function(x) length(unique(x)))
names(gene_rank)[[2]] <- "Distinct_pathways"; gene_rank <- gene_rank[order(-gene_rank$Distinct_pathways, gene_rank$Gene), ]
keep_genes <- head(gene_rank$Gene, max_genes)
display <- triples[triples$Ingredient %in% keep_ingredients & triples$Gene %in% keep_genes, ]
write_out(display, "02_sankey_display_triples.csv")
node_meta <- unique(rbind(
  data.frame(term = display$Ingredient, type = "Ingredient", label = display$Ingredient_name,
    pvalue = NA_real_, minus_log10_p = NA_real_),
  data.frame(term = display$Gene, type = "Gene", label = display$Gene,
    pvalue = NA_real_, minus_log10_p = NA_real_),
  data.frame(term = display$Pathway, type = "Pathway", label = display$Pathway_name,
    pvalue = display$Pathway_pvalue, minus_log10_p = display$Pathway_minus_log10_p)
))
write_out(node_meta, "03_sankey_nodes.csv")

plot_status <- "not attempted"
if (requireNamespace("ggsankeyfier", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggsankeyfier))
  wide <- display[, c("Ingredient", "Gene", "Pathway", "weight")]
  pathway_label <- unique(display[, c("Pathway", "Pathway_name")])
  pathway_lookup <- setNames(pathway_label$Pathway_name, pathway_label$Pathway)
  labels <- unname(pathway_lookup[wide$Pathway])
  labels[is.na(labels) | labels == ""] <- wide$Pathway[is.na(labels) | labels == ""]
  wide$Pathway <- ifelse(nchar(labels) > 44, paste0(substr(labels, 1, 43), "..."), labels)
  long <- ggsankeyfier::pivot_stages_longer(wide, stages_from = c("Ingredient", "Gene", "Pathway"), values_from = "weight")
  p <- ggplot2::ggplot(long, ggplot2::aes(x = stage, y = weight, group = node, edge_id = edge_id, connector = connector)) +
    ggsankeyfier::geom_sankeyedge(ggplot2::aes(fill = stage), alpha = 0.45,
      position = ggsankeyfier::position_sankey(v_space = "auto", width = 0.16)) +
    ggsankeyfier::geom_sankeynode(ggplot2::aes(fill = stage), color = "grey25",
      position = ggsankeyfier::position_sankey(v_space = "auto", width = 0.16)) +
    ggplot2::geom_text(ggplot2::aes(label = node), stat = "sankeynode", size = 2.4,
      position = ggsankeyfier::position_sankey(v_space = "auto", width = 0.16), check_overlap = TRUE) +
    ggplot2::scale_fill_manual(values = c(Ingredient = "#D9A441", Gene = "#4C78A8", Pathway = "#59A14F")) +
    ggplot2::labs(x = NULL, y = "Unique ingredient-gene-pathway relations") +
    ggplot2::theme_void(base_size = 10) + ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(10, 20, 10, 20))
  height <- max(8, min(20, 0.16 * length(unique(c(display$Ingredient, display$Gene, display$Pathway))) + 4))
  ggplot2::ggsave(file.path(output_dir, "04_ingredient_gene_pathway_sankey.png"), p, width = 11, height = height, dpi = 600, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "04_ingredient_gene_pathway_sankey.pdf"), p, width = 11, height = height, bg = "white")
  plot_status <- "generated PNG 600 dpi and vector PDF"
} else {
  warning("Sankey tables were generated, but plotting requires CRAN packages ggsankeyfier and ggplot2.")
  plot_status <- "tables only; ggsankeyfier or ggplot2 unavailable"
}
readability_warning <- if (nrow(display) > 2000) "Dense figure: retain full tables and consider a separately labeled smaller display-only rerun" else "none"
qc <- data.frame(Item = c("All triples", "Displayed triples", "Displayed ingredients", "Displayed genes", "Displayed pathways", "Plot status", "Readability warning"),
  Value = c(nrow(triples), nrow(display), length(unique(display$Ingredient)), length(unique(display$Gene)), length(unique(display$Pathway)), plot_status, readability_warning))
write_out(qc, "05_sankey_qc.csv")
cat("all_triples=", nrow(triples), "\n", sep = "")
cat("display_triples=", nrow(display), "\n", sep = "")
cat("plot_status=", plot_status, "\n", sep = "")
