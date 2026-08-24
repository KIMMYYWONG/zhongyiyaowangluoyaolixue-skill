if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop(paste(
    "Usage: Rscript 02b_disease_comorbidity_intersection.R",
    "<disease_a_targets.csv|tsv|txt> <disease_b_targets.csv|tsv|txt>",
    "<disease_a_tracking.csv|tsv|txt|NONE> <disease_b_tracking.csv|tsv|txt|NONE>",
    "<output_dir> <disease_a_label> <disease_b_label> [set_label=broad_candidate]"
  ))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

a_file <- args[[1]]
b_file <- args[[2]]
a_tracking_file <- args[[3]]
b_tracking_file <- args[[4]]
output_dir <- args[[5]]
a_label <- trimws(args[[6]])
b_label <- trimws(args[[7]])
set_label <- if (length(args) >= 8) trimws(args[[8]]) else "broad_candidate"

if (a_label == "" || b_label == "" || a_label == b_label) {
  stop("Disease labels must be non-empty and distinct.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  signature <- readBin(path, what = "raw", n = 4)
  if (length(signature) >= 2 && identical(as.integer(signature[1:2]), c(80L, 75L))) {
    stop("Input is an XLSX/ZIP file but has a delimited-table workflow path; export it as CSV/TSV or use an explicit XLSX adapter: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (ext %in% c("tsv", "txt")) return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  stop("Only CSV, TSV, and TXT inputs are supported: ", path)
}

normalize_key <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
find_gene_column <- function(df) {
  keys <- normalize_key(names(df))
  wanted <- normalize_key(c("Gene_symbol", "Gene symbol", "Gene.Symbol", "symbol", "gene", "target"))
  hit <- match(wanted, keys, nomatch = 0)
  hit <- hit[hit > 0]
  if (!length(hit)) stop("No gene-symbol column found. Available columns: ", paste(names(df), collapse = ", "))
  names(df)[hit[[1]]]
}
clean_genes <- function(df) {
  col <- find_gene_column(df)
  genes <- toupper(trimws(as.character(df[[col]])))
  genes <- genes[!is.na(genes) & genes != "" & grepl("^[A-Z0-9][A-Z0-9._-]*$", genes)]
  sort(unique(genes))
}

a_raw <- read_table(a_file)
b_raw <- read_table(b_file)
a_genes <- clean_genes(a_raw)
b_genes <- clean_genes(b_raw)
if (!length(a_genes) || !length(b_genes)) stop("Each disease target set must contain at least one valid gene symbol.")

union_genes <- sort(union(a_genes, b_genes))
intersection_genes <- sort(intersect(a_genes, b_genes))
if (!length(intersection_genes)) stop("The two disease target sets have an empty intersection; review disease identity, set choice, and identifiers.")

membership <- data.frame(
  Gene_symbol = union_genes,
  disease_a = union_genes %in% a_genes,
  disease_b = union_genes %in% b_genes,
  membership = ifelse(union_genes %in% intersection_genes, "both", ifelse(union_genes %in% a_genes, "disease_a_only", "disease_b_only")),
  disease_a_label = a_label,
  disease_b_label = b_label,
  input_set = set_label,
  stringsAsFactors = FALSE
)

intersection_table <- data.frame(
  Gene_symbol = intersection_genes,
  Disease_A = a_label,
  Disease_B = b_label,
  Input_set = set_label,
  stringsAsFactors = FALSE
)

read_tracking <- function(path, disease_label) {
  if (toupper(path) == "NONE") return(data.frame())
  x <- read_table(path)
  gene_col <- find_gene_column(x)
  x$Gene_symbol <- toupper(trimws(as.character(x[[gene_col]])))
  x <- x[x$Gene_symbol %in% intersection_genes, , drop = FALSE]
  x$Disease_label <- disease_label
  x$Input_set <- set_label
  x
}

tracking <- bind_rows(read_tracking(a_tracking_file, a_label), read_tracking(b_tracking_file, b_label))
if (!nrow(tracking)) {
  tracking <- intersection_table %>% transmute(Gene_symbol, Disease_label = paste(a_label, b_label, sep = "; "), Input_set)
}

stats <- data.frame(
  Measure = c("Disease A unique targets", "Disease B unique targets", "Union targets", "Comorbidity intersection targets", "Jaccard index"),
  Value = c(length(a_genes), length(b_genes), length(union_genes), length(intersection_genes), length(intersection_genes) / length(union_genes)),
  Disease_A = a_label,
  Disease_B = b_label,
  Input_set = set_label,
  stringsAsFactors = FALSE
)

qc <- data.frame(
  Check = c("Distinct non-empty disease labels", "Disease A target set nonempty", "Disease B target set nonempty", "Comorbidity intersection nonempty", "Intersection is subset of disease A", "Intersection is subset of disease B", "Membership gene keys unique"),
  Pass = c(a_label != "" && b_label != "" && a_label != b_label, length(a_genes) > 0, length(b_genes) > 0, length(intersection_genes) > 0, all(intersection_genes %in% a_genes), all(intersection_genes %in% b_genes), !anyDuplicated(membership$Gene_symbol)),
  Detail = c(paste(a_label, b_label, sep = " vs "), length(a_genes), length(b_genes), length(intersection_genes), "", "", length(union_genes)),
  stringsAsFactors = FALSE
)

write.csv(membership, file.path(output_dir, "01_疾病集合成员表.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(intersection_table, file.path(output_dir, "02_共病交集靶点.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tracking, file.path(output_dir, "03_共病交集来源追踪表.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stats, file.path(output_dir, "04_共病交集统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(qc, file.path(output_dir, "05_共病交集QC.csv"), row.names = FALSE, fileEncoding = "UTF-8")

plot_df <- data.frame(
  Category = factor(c(paste0(a_label, " only"), "Intersection", paste0(b_label, " only")), levels = c(paste0(a_label, " only"), "Intersection", paste0(b_label, " only"))),
  Count = c(length(setdiff(a_genes, b_genes)), length(intersection_genes), length(setdiff(b_genes, a_genes)))
)
p <- ggplot(plot_df, aes(Category, Count, fill = Category)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.35, size = 4) +
  scale_fill_manual(values = c("#4C78A8", "#7A5195", "#F58518")) +
  labs(title = paste(a_label, "and", b_label, "target overlap"), subtitle = paste("Input set:", set_label), x = NULL, y = "Unique gene count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(output_dir, "06_共病靶点重叠图.png"), p, width = 8, height = 5.5, dpi = 600)
ggsave(file.path(output_dir, "06_共病靶点重叠图.pdf"), p, width = 8, height = 5.5)

if (!all(qc$Pass)) stop("Comorbidity intersection QC failed; inspect 05_共病交集QC.csv.")
cat("disease_a_targets=", length(a_genes), "\n", sep = "")
cat("disease_b_targets=", length(b_genes), "\n", sep = "")
cat("comorbidity_intersection_targets=", length(intersection_genes), "\n", sep = "")

