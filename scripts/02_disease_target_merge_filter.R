if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(paste(
    "Usage: Rscript 02_disease_target_merge_filter.R <project_dir>",
    "[raw_dir] [output_dir] [policy_csv]"
  ))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
raw_dir <- if (length(args) >= 2) args[[2]] else file.path(project_dir, "疾病靶点预测")
output_dir <- if (length(args) >= 3) args[[3]] else file.path(project_dir, "01.Result", "02.疾病靶点汇总")
policy_file <- if (length(args) >= 4) args[[4]] else file.path(project_dir, "disease-source-policy.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_data_file <- function(path) {
  signature <- readBin(path, what = "raw", n = 4)
  if (length(signature) >= 2 && identical(as.integer(signature[1:2]), c(80L, 75L))) {
    stop("Input is an XLSX/ZIP file with a delimited-table extension; rename or export it correctly before screening: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  }
  if (ext %in% c("tsv", "txt")) {
    return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  }
  stop("Only CSV, TSV, and TXT inputs are supported: ", path)
}

normalize_key <- function(x) {
  gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
}

find_column <- function(df, candidates) {
  keys <- normalize_key(names(df))
  candidate_keys <- normalize_key(candidates)
  hit <- match(candidate_keys, keys, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) == 0) return(NA_character_)
  names(df)[hit[[1]]]
}

collapse_text_columns <- function(df, columns) {
  columns <- unique(columns[!is.na(columns) & columns %in% names(df)])
  if (length(columns) == 0) return(rep("", nrow(df)))
  apply(df[, columns, drop = FALSE], 1, function(row) {
    values <- trimws(as.character(row))
    values <- unique(values[!is.na(values) & values != "" & values != "NA"])
    paste(values, collapse = " | ")
  })
}

parse_logical <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value %in% c("true", "t", "1", "yes", "y", "curated", "direct")
}

collapse_unique <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != "" & x != "NA"]
  paste(sort(x), collapse = "; ")
}

if (!dir.exists(raw_dir)) stop("Disease target input directory not found: ", raw_dir)
input_files <- list.files(raw_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE, ignore.case = TRUE)
if (length(input_files) == 0) stop("No CSV, TSV, or TXT disease target files found in: ", raw_dir)
if (!file.exists(policy_file)) stop("Disease source policy file not found: ", policy_file)

policy <- read.csv(policy_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
required_policy_columns <- c("Data_base", "screen_mode", "quantile", "min_score", "curated_keywords")
if (!all(required_policy_columns %in% names(policy))) {
  stop("Policy file must contain: ", paste(required_policy_columns, collapse = ", "))
}
policy$Source_key <- normalize_key(policy$Data_base)

allowed_modes <- c("median_score", "curated_or_median_score", "curated_all", "all", "min_score")
if (any(policy$Source_key == "") || anyDuplicated(policy$Source_key)) {
  stop("Disease source policy contains empty or duplicate Data_base values after normalization.")
}
invalid_modes <- unique(policy$screen_mode[!policy$screen_mode %in% allowed_modes])
if (length(invalid_modes)) {
  stop("Unsupported disease source screen_mode: ", paste(invalid_modes, collapse = ", "))
}
median_rows <- policy$screen_mode %in% c("median_score", "curated_or_median_score")
median_quantiles <- suppressWarnings(as.numeric(policy$quantile[median_rows]))
if (length(median_quantiles) && any(!is.finite(median_quantiles) | median_quantiles < 0 | median_quantiles > 1)) {
  stop("Median-based disease source policies require quantile values in [0, 1].")
}
minimum_rows <- policy$screen_mode == "min_score"
minimum_scores <- suppressWarnings(as.numeric(policy$min_score[minimum_rows]))
if (length(minimum_scores) && any(!is.finite(minimum_scores))) {
  stop("min_score disease source policies require a finite min_score value.")
}

input_source_keys <- normalize_key(tools::file_path_sans_ext(basename(input_files)))
missing_policy_sources <- unique(tools::file_path_sans_ext(basename(input_files))[!input_source_keys %in% policy$Source_key])
if (length(missing_policy_sources)) {
  stop("Missing disease source policy for input file source(s): ", paste(missing_policy_sources, collapse = ", "))
}

all_records <- list()
for (path in input_files) {
  df <- read_data_file(path)
  if (nrow(df) == 0) next

  database_name <- tools::file_path_sans_ext(basename(path))
  gene_col <- find_column(df, c(
    "Gene_symbol", "Gene symbol", "Symbol", "Approved_symbol",
    "Approved Symbol", "Target_symbol", "target", "gene"
  ))
  if (is.na(gene_col)) stop("Missing gene-symbol column in: ", path)

  score_col <- find_column(df, c(
    "Score", "Relevance_score", "Relevance score", "Inference_score",
    "Inference score", "Overall_association_score", "overall score",
    "GDA_score", "DISGENET_score"
  ))
  curated_col <- find_column(df, c("Is_curated", "Curated", "Direct_evidence", "Is_direct"))
  evidence_cols <- c(
    find_column(df, c("Evidence_type", "Evidence type")),
    find_column(df, c("Association_type", "Association type")),
    find_column(df, c("Evidence", "Relation", "Relationship", "Action")),
    find_column(df, c("Source", "Evidence_source", "Evidence source"))
  )

  gene_raw <- trimws(as.character(df[[gene_col]]))
  gene_symbol <- toupper(gene_raw)
  gene_symbol[gene_symbol %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  score_numeric <- if (is.na(score_col)) rep(NA_real_, nrow(df)) else suppressWarnings(as.numeric(df[[score_col]]))
  explicit_curated <- if (is.na(curated_col)) rep(FALSE, nrow(df)) else parse_logical(df[[curated_col]])
  evidence_text <- collapse_text_columns(df, evidence_cols)

  raw_copy <- df
  names(raw_copy) <- paste0("raw__", make.unique(names(raw_copy)))
  canonical <- data.frame(
    Record_id = sprintf("%s:%06d", database_name, seq_len(nrow(df))),
    Data_base = database_name,
    Source_file = basename(path),
    Source_key = normalize_key(database_name),
    Gene_symbol_raw = gene_raw,
    Gene_symbol = gene_symbol,
    Score_field = if (is.na(score_col)) "" else score_col,
    Score_numeric = score_numeric,
    Evidence_text = evidence_text,
    Explicit_curated = explicit_curated,
    stringsAsFactors = FALSE
  )
  all_records[[length(all_records) + 1]] <- bind_cols(canonical, raw_copy)
}

if (length(all_records) == 0) stop("All disease target input files were empty.")
combined <- bind_rows(all_records)
combined$Valid_gene <- !is.na(combined$Gene_symbol) & grepl("^[A-Z0-9][A-Z0-9._-]*$", combined$Gene_symbol)
combined$Screen_mode <- ""
combined$Source_cutoff <- NA_real_
combined$Evidence_preferred <- FALSE
combined$Source_pass <- FALSE
combined$Screen_reason <- ""

source_keys <- unique(combined$Source_key)
for (source_key in source_keys) {
  idx <- which(combined$Source_key == source_key)
  policy_idx <- which(policy$Source_key == source_key)
  if (length(policy_idx) != 1) {
    combined$Screen_reason[idx] <- "missing or ambiguous source policy"
    next
  }

  rule <- policy[policy_idx, , drop = FALSE]
  mode <- trimws(as.character(rule$screen_mode[[1]]))
  q <- suppressWarnings(as.numeric(rule$quantile[[1]]))
  min_score <- suppressWarnings(as.numeric(rule$min_score[[1]]))
  keywords <- trimws(as.character(rule$curated_keywords[[1]]))
  keyword_match <- if (is.na(keywords) || keywords == "") {
    rep(FALSE, length(idx))
  } else {
    grepl(keywords, combined$Evidence_text[idx], ignore.case = TRUE, perl = TRUE)
  }
  preferred <- combined$Explicit_curated[idx] | keyword_match
  if (mode == "curated_all") preferred <- combined$Valid_gene[idx]

  scores <- combined$Score_numeric[idx]
  finite_scores <- scores[is.finite(scores)]
  cutoff <- NA_real_
  if (mode %in% c("median_score", "curated_or_median_score")) {
    if (!is.finite(q)) q <- 0.5
    if (length(finite_scores) > 0) cutoff <- as.numeric(quantile(finite_scores, probs = q, na.rm = TRUE, names = FALSE, type = 7))
  } else if (mode == "min_score") {
    cutoff <- min_score
  }

  score_pass <- is.finite(scores) & is.finite(cutoff) & scores >= cutoff
  pass <- switch(
    mode,
    median_score = score_pass,
    curated_or_median_score = preferred | score_pass,
    curated_all = combined$Valid_gene[idx],
    all = combined$Valid_gene[idx],
    min_score = score_pass,
    rep(FALSE, length(idx))
  )
  pass <- pass & combined$Valid_gene[idx]

  reason <- ifelse(
    !combined$Valid_gene[idx], "invalid or missing standardized gene symbol",
    ifelse(
      pass & preferred, "retained by curated or direct evidence",
      ifelse(
        pass & score_pass, sprintf("retained by %s cutoff %.6g", mode, cutoff),
        ifelse(
          mode %in% c("median_score", "curated_or_median_score", "min_score") & !is.finite(scores),
          "excluded: required score missing",
          ifelse(mode %in% c("median_score", "curated_or_median_score", "min_score"),
                 sprintf("excluded: score below cutoff %.6g", cutoff),
                 ifelse(pass, "retained by source policy", "excluded by source policy"))
        )
      )
    )
  )

  combined$Screen_mode[idx] <- mode
  combined$Source_cutoff[idx] <- cutoff
  combined$Evidence_preferred[idx] <- preferred
  combined$Source_pass[idx] <- pass
  combined$Screen_reason[idx] <- reason
}

raw_output <- combined %>%
  select(-Screen_mode, -Source_cutoff, -Evidence_preferred, -Source_pass, -Screen_reason)
decision_output <- combined

source_summary <- decision_output %>%
  group_by(Data_base, Source_key, Screen_mode) %>%
  summarise(
    Raw_records = n(),
    Raw_unique_targets = n_distinct(Gene_symbol[Valid_gene]),
    Scored_records = sum(is.finite(Score_numeric)),
    Source_cutoff = if (all(is.na(Source_cutoff))) NA_real_ else first(na.omit(Source_cutoff)),
    Preferred_evidence_records = sum(Evidence_preferred & Valid_gene),
    Retained_records = sum(Source_pass),
    Retained_unique_targets = n_distinct(Gene_symbol[Source_pass]),
    Missing_gene_records = sum(!Valid_gene),
    .groups = "drop"
  )

tracking <- decision_output %>%
  filter(Source_pass) %>%
  group_by(Gene_symbol, Data_base) %>%
  summarise(
    Source_record_count = n(),
    Score_max = if (all(!is.finite(Score_numeric))) NA_real_ else max(Score_numeric, na.rm = TRUE),
    Source_cutoff = if (all(is.na(Source_cutoff))) NA_real_ else first(na.omit(Source_cutoff)),
    Evidence_preferred = any(Evidence_preferred),
    Evidence_summary = collapse_unique(Evidence_text),
    Record_ids = collapse_unique(Record_id),
    .groups = "drop"
  )

gene_summary <- tracking %>%
  group_by(Gene_symbol) %>%
  summarise(
    Data_bases = collapse_unique(Data_base),
    Database_count = n_distinct(Data_base),
    Curated_or_direct = any(Evidence_preferred),
    Source_record_count = sum(Source_record_count),
    Score_max_across_sources = if (all(!is.finite(Score_max))) NA_real_ else max(Score_max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Evidence_tier = case_when(
      Curated_or_direct ~ "A_curated_or_direct",
      Database_count >= 2 ~ "B_multi_database",
      TRUE ~ "C_single_source_screened"
    ),
    Broad_candidate = TRUE,
    Recommended_include = Curated_or_direct | Database_count >= 2
  ) %>%
  arrange(Evidence_tier, desc(Database_count), Gene_symbol)

broad_targets <- gene_summary %>% filter(Broad_candidate)
recommended_targets <- gene_summary %>% filter(Recommended_include)

if (nrow(broad_targets) == 0) {
  stop("Disease target screening retained no valid targets; review source policies and input identifiers.")
}

write.csv(raw_output, file.path(output_dir, "01_疾病靶点原始汇总.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decision_output, file.path(output_dir, "02_疾病靶点筛选决策.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tracking, file.path(output_dir, "03_疾病靶点来源追踪.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(broad_targets, file.path(output_dir, "04_疾病靶点候选全集.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(recommended_targets, file.path(output_dir, "05_疾病靶点推荐集.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(source_summary, file.path(output_dir, "06_疾病靶点筛选统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")

plot_df <- bind_rows(
  source_summary %>% transmute(Data_base, Target_set = "Raw", Count = Raw_unique_targets),
  source_summary %>% transmute(Data_base, Target_set = "Retained", Count = Retained_unique_targets)
)
p <- ggplot(plot_df, aes(x = Data_base, y = Count, fill = Target_set)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.75), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(Raw = "#94A3B8", Retained = "#2563EB")) +
  labs(x = NULL, y = "Unique targets", fill = "Target set") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"))
ggsave(file.path(output_dir, "07_疾病靶点数据库统计.png"), p, width = 9, height = 6, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "07_疾病靶点数据库统计.pdf"), p, width = 9, height = 6, bg = "white")

cat("raw_records=", nrow(raw_output), "\n", sep = "")
cat("broad_unique_targets=", nrow(broad_targets), "\n", sep = "")
cat("recommended_unique_targets=", nrow(recommended_targets), "\n", sep = "")
cat("output_dir=", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
