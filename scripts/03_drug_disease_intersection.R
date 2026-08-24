if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: Rscript 03_drug_disease_intersection.R <project_dir>",
    "<drug_target_file> <disease_target_file> <disease_tracking_file>",
    "[output_dir] [drug_label] [disease_label] [disease_set_label]"
  ))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(VennDiagram)
  library(grid)
})

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
resolve_path <- function(path) {
  if (grepl("^([A-Za-z]:|/)", path)) path else file.path(project_dir, path)
}

drug_file <- resolve_path(args[[2]])
disease_file <- resolve_path(args[[3]])
disease_tracking_file <- resolve_path(args[[4]])
output_dir <- if (length(args) >= 5) resolve_path(args[[5]]) else file.path(project_dir, "01.Result", "03.药物-疾病 交集靶点")
drug_label <- if (length(args) >= 6) args[[6]] else "Drug targets"
disease_label <- if (length(args) >= 7) args[[7]] else "Disease targets"
disease_set_label <- if (length(args) >= 8) args[[8]] else "recommended"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_data_file <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (ext %in% c("tsv", "txt")) return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  stop("Only CSV, TSV, and TXT inputs are supported: ", path)
}

normalize_key <- function(x) {
  gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
}

find_column <- function(df, candidates, required = FALSE) {
  keys <- normalize_key(names(df))
  candidate_keys <- normalize_key(candidates)
  hit <- match(candidate_keys, keys, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) == 0) {
    if (required) stop("Missing required column. Expected one of: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  names(df)[hit[[1]]]
}

collapse_unique <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != "" & x != "NA"]
  paste(sort(x), collapse = "; ")
}

prepare_gene_column <- function(df) {
  gene_col <- find_column(df, c(
    "Gene_symbol", "Gene symbol", "gene_symbol_primary", "Symbol",
    "Approved_symbol", "target"
  ), required = TRUE)
  df$Gene_symbol <- toupper(trimws(as.character(df[[gene_col]])))
  df <- df[!is.na(df$Gene_symbol) & df$Gene_symbol != "" & df$Gene_symbol != "NA", , drop = FALSE]
  df
}

drug <- prepare_gene_column(read_data_file(drug_file))
disease <- prepare_gene_column(read_data_file(disease_file))
disease_tracking <- prepare_gene_column(read_data_file(disease_tracking_file))

drug_source_col <- find_column(drug, c("source", "Data_base", "Drug_source", "database"))
herb_col <- find_column(drug, c("herb_normalized_cn", "herb_cn", "Herb", "herb_pinyin"))
compound_col <- find_column(drug, c("compound_name_preferred", "Molecule_name", "Compound", "Mol_ID", "tcmsp_molid"))

drug$Drug_source_value <- if (is.na(drug_source_col)) "unspecified" else as.character(drug[[drug_source_col]])
drug$Herb_value <- if (is.na(herb_col)) "" else as.character(drug[[herb_col]])
drug$Compound_value <- if (is.na(compound_col)) "" else as.character(drug[[compound_col]])

drug_summary <- drug %>%
  group_by(Gene_symbol) %>%
  summarise(
    Drug_sources = collapse_unique(Drug_source_value),
    Herbs = collapse_unique(Herb_value),
    Compounds = collapse_unique(Compound_value),
    Drug_record_count = n(),
    .groups = "drop"
  )

disease_source_col <- find_column(disease_tracking, c("Data_base", "source", "database"))
evidence_col <- find_column(disease_tracking, c("Evidence_summary", "Evidence_text", "Evidence_type"))
score_col <- find_column(disease_tracking, c("Score_max", "Score_numeric", "Score"))
disease_tracking$Disease_source_value <- if (is.na(disease_source_col)) "unspecified" else as.character(disease_tracking[[disease_source_col]])
disease_tracking$Evidence_value <- if (is.na(evidence_col)) "" else as.character(disease_tracking[[evidence_col]])
disease_tracking$Score_value <- if (is.na(score_col)) NA_real_ else suppressWarnings(as.numeric(disease_tracking[[score_col]]))

disease_summary <- disease_tracking %>%
  filter(Gene_symbol %in% disease$Gene_symbol) %>%
  group_by(Gene_symbol) %>%
  summarise(
    Disease_databases = collapse_unique(Disease_source_value),
    Disease_database_count = n_distinct(Disease_source_value),
    Disease_evidence = collapse_unique(Evidence_value),
    Disease_score_max = if (all(!is.finite(Score_value))) NA_real_ else max(Score_value, na.rm = TRUE),
    Disease_record_count = n(),
    .groups = "drop"
  )

disease_meta_cols <- intersect(c("Gene_symbol", "Data_bases", "Database_count", "Curated_or_direct", "Evidence_tier", "Recommended_include"), names(disease))
disease_meta <- disease[, disease_meta_cols, drop = FALSE] %>% distinct(Gene_symbol, .keep_all = TRUE)
disease_summary <- full_join(disease_summary, disease_meta, by = "Gene_symbol")

drug_targets <- sort(unique(drug_summary$Gene_symbol))
disease_targets <- sort(unique(disease$Gene_symbol))
all_targets <- sort(union(drug_targets, disease_targets))
intersection_targets <- sort(intersect(drug_targets, disease_targets))

if (length(intersection_targets) == 0) {
  stop(
    "Drug-disease intersection is empty (drug targets=", length(drug_targets),
    ", disease targets=", length(disease_targets),
    "). Review the frozen target sets before PPI or enrichment."
  )
}

membership <- data.frame(
  Gene_symbol = all_targets,
  In_drug_set = all_targets %in% drug_targets,
  In_disease_set = all_targets %in% disease_targets,
  Set_membership = ifelse(
    all_targets %in% intersection_targets,
    "intersection",
    ifelse(all_targets %in% drug_targets, "drug_only", "disease_only")
  ),
  stringsAsFactors = FALSE
) %>%
  left_join(drug_summary, by = "Gene_symbol") %>%
  left_join(disease_summary, by = "Gene_symbol")

intersection_table <- data.frame(Gene_symbol = intersection_targets, stringsAsFactors = FALSE)
intersection_trace <- membership %>% filter(Set_membership == "intersection")

write.csv(membership, file.path(output_dir, "01_药物疾病集合成员表.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(intersection_table, file.path(output_dir, "02_药物疾病交集靶点.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(intersection_trace, file.path(output_dir, "03_交集靶点来源追踪表.csv"), row.names = FALSE, fileEncoding = "UTF-8")

venn_values <- list()
venn_values[[drug_label]] <- drug_targets
venn_values[[disease_label]] <- disease_targets

venn_png <- file.path(output_dir, "04_药物疾病靶点韦恩图.png")
venn_pdf <- file.path(output_dir, "04_药物疾病靶点韦恩图.pdf")
venn.diagram(
  x = venn_values,
  filename = venn_png,
  imagetype = "png",
  height = 3000,
  width = 3200,
  resolution = 500,
  units = "px",
  fill = c("#FBBF24", "#7C3AED"),
  alpha = c(0.58, 0.52),
  col = NA,
  cat.col = c("#111827", "#111827"),
  cat.cex = 1.2,
  cex = 1.25,
  scaled = FALSE
)

if (capabilities("cairo")) {
  cairo_pdf(venn_pdf, width = 8, height = 7, family = "Microsoft YaHei")
} else {
  pdf(venn_pdf, width = 8, height = 7)
}
grid.draw(venn.diagram(
  x = venn_values,
  filename = NULL,
  fill = c("#FBBF24", "#7C3AED"),
  alpha = c(0.58, 0.52),
  col = NA,
  cat.col = c("#111827", "#111827"),
  cat.fontfamily = "Microsoft YaHei",
  cat.cex = 1.2,
  fontfamily = "Microsoft YaHei",
  cex = 1.25,
  scaled = FALSE
))
dev.off()

statistics <- data.frame(
  Metric = c("drug_unique_targets", "disease_unique_targets", "union_targets", "intersection_targets", "disease_set_label"),
  Value = c(length(drug_targets), length(disease_targets), length(all_targets), length(intersection_targets), disease_set_label),
  stringsAsFactors = FALSE
)
write.csv(statistics, file.path(output_dir, "05_交集统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("drug_unique_targets=", length(drug_targets), "\n", sep = "")
cat("disease_unique_targets=", length(disease_targets), "\n", sep = "")
cat("intersection_targets=", length(intersection_targets), "\n", sep = "")
cat("output_dir=", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
