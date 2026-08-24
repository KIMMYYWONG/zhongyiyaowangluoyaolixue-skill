args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript prepare_legacy_two_herbs.R <herb_targets.xlsx> <output_dir>")

suppressPackageStartupMessages(library(openxlsx))
input_file <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clean <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | toupper(z) %in% c("N/A", "NA", "NULL")] <- ""
  z
}

read_legacy_sheet <- function(sheet_index, herb_cn, herb_pinyin) {
  x <- read.xlsx(input_file, sheet = sheet_index, colNames = FALSE, skipEmptyRows = TRUE)
  if (ncol(x) < 7L) stop("Legacy sheet has fewer than seven columns: ", sheet_index)
  x <- x[, 1:7]
  names(x) <- c("herb_pinyin_raw", "pubchem_cid_raw", "source_ingredient_id",
                "compound_name_raw", "gene_symbol_input", "cas_raw", "source")
  x[] <- lapply(x, clean)
  x$herb_normalized_cn <- herb_cn
  x$herb_pinyin <- herb_pinyin
  x$source <- toupper(x$source)
  x$source[x$source == "HERB"] <- "HERB2"
  x$legacy_row <- seq_len(nrow(x))
  x <- x[x$source_ingredient_id != "" & x$compound_name_raw != "" & x$gene_symbol_input != "", ]
  x[, c("herb_normalized_cn", "herb_pinyin", "legacy_row", "source",
        "source_ingredient_id", "compound_name_raw", "pubchem_cid_raw",
        "cas_raw", "gene_symbol_input")]
}

# Workbook sheet order was inspected before extraction: Bai Shao is sheet 1 and Xi Xin is sheet 15.
# Construct Chinese labels from code points so the script remains reproducible under Windows R locales.
bai_shao_cn <- intToUtf8(c(0x767D, 0x828D))
xi_xin_cn <- intToUtf8(c(0x7EC6, 0x8F9B))
bai_shao <- read_legacy_sheet(1, bai_shao_cn, "Bai Shao")
xi_xin <- read_legacy_sheet(15, xi_xin_cn, "Xi Xin")
raw <- rbind(bai_shao, xi_xin)

relation_cols <- c("herb_normalized_cn", "source", "source_ingredient_id",
                   "compound_name_raw", "pubchem_cid_raw", "cas_raw", "gene_symbol_input")
raw$exact_duplicate <- duplicated(raw[, relation_cols])
dedup <- raw[!raw$exact_duplicate, ]

write.csv(raw, file.path(output_dir, "00_legacy_two_herb_rows.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(dedup, file.path(output_dir, "01_legacy_exact_deduplicated.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summarize_herb <- function(df) {
  data.frame(
    herb_normalized_cn = df$herb_normalized_cn[[1]],
    raw_relationship_rows = nrow(df),
    exact_duplicate_rows = sum(df$exact_duplicate),
    unique_source_ingredient_ids = length(unique(df$source_ingredient_id)),
    unique_pubchem_cids = length(unique(df$pubchem_cid_raw[df$pubchem_cid_raw != ""])),
    rows_missing_pubchem_cid = sum(df$pubchem_cid_raw == ""),
    rows_missing_cas = sum(df$cas_raw == ""),
    unique_input_gene_symbols = length(unique(df$gene_symbol_input)),
    sources = paste(sort(unique(df$source)), collapse = "|"),
    stringsAsFactors = FALSE
  )
}
summary <- do.call(rbind, lapply(split(raw, raw$herb_pinyin), summarize_herb))
write.csv(summary, file.path(output_dir, "02_preflight_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

compound_key_audit <- aggregate(
  cbind(pubchem_count = raw$pubchem_cid_raw, name_count = raw$compound_name_raw, cas_count = raw$cas_raw),
  by = raw[c("herb_normalized_cn", "source", "source_ingredient_id")],
  FUN = function(x) length(unique(x[x != ""]))
)
compound_key_audit$identifier_conflict <- compound_key_audit$pubchem_count > 1L |
  compound_key_audit$name_count > 1L | compound_key_audit$cas_count > 1L
write.csv(compound_key_audit, file.path(output_dir, "03_source_ingredient_identity_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

gene_symbols <- data.frame(target_name_raw = sort(unique(toupper(dedup$gene_symbol_input))), stringsAsFactors = FALSE)
write.csv(gene_symbols, file.path(output_dir, "04_unique_gene_symbols.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("raw_rows=", nrow(raw), "\n", sep = "")
cat("deduplicated_rows=", nrow(dedup), "\n", sep = "")
cat("source_ingredient_conflicts=", sum(compound_key_audit$identifier_conflict), "\n", sep = "")
