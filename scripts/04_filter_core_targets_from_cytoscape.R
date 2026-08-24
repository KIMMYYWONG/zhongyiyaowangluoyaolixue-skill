if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste(
    "Usage: Rscript 04_filter_core_targets_from_cytoscape.R",
    "<cytoscape_node_table.csv|tsv> <output_dir> [node_name_column]"
  ))
}

options(stringsAsFactors = FALSE)

input_file <- args[[1]]
output_dir <- args[[2]]
requested_node_col <- if (length(args) >= 3) args[[3]] else ""
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_data_file <- function(path) {
  if (!file.exists(path)) stop("Cytoscape node table not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (ext %in% c("tsv", "txt")) return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
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

nodes <- read_data_file(input_file)
if (nrow(nodes) == 0) stop("Cytoscape node table is empty.")

node_col <- if (requested_node_col != "" && requested_node_col %in% names(nodes)) {
  requested_node_col
} else {
  find_column(nodes, c("Gene_symbol", "name", "shared name", "display name", "query term", "preferredName"))
}
betweenness_col <- find_column(nodes, c("BetweennessCentrality", "Betweenness Centrality", "Betweenness"))
closeness_col <- find_column(nodes, c("ClosenessCentrality", "Closeness Centrality", "Closeness"))
degree_col <- find_column(nodes, c("Degree", "degree.layout"))

missing_columns <- c(
  if (is.na(node_col)) "node name" else character(),
  if (is.na(betweenness_col)) "Betweenness" else character(),
  if (is.na(closeness_col)) "Closeness" else character(),
  if (is.na(degree_col)) "Degree" else character()
)
if (length(missing_columns) > 0) {
  stop("Missing required Cytoscape node-table columns: ", paste(missing_columns, collapse = ", "))
}

decision <- nodes
decision$Gene_symbol_raw <- trimws(as.character(nodes[[node_col]]))
decision$Gene_symbol <- toupper(decision$Gene_symbol_raw)
decision$Betweenness_value <- suppressWarnings(as.numeric(nodes[[betweenness_col]]))
decision$Closeness_value <- suppressWarnings(as.numeric(nodes[[closeness_col]]))
decision$Degree_value <- suppressWarnings(as.numeric(nodes[[degree_col]]))
decision$Valid_metrics <- is.finite(decision$Betweenness_value) &
  is.finite(decision$Closeness_value) &
  is.finite(decision$Degree_value) &
  !is.na(decision$Gene_symbol) & decision$Gene_symbol != ""

if (!any(decision$Valid_metrics)) stop("No Cytoscape nodes have all three valid topology metrics.")

betweenness_median <- median(decision$Betweenness_value[decision$Valid_metrics], na.rm = TRUE)
closeness_median <- median(decision$Closeness_value[decision$Valid_metrics], na.rm = TRUE)
degree_median <- median(decision$Degree_value[decision$Valid_metrics], na.rm = TRUE)

decision$Betweenness_above_median <- decision$Valid_metrics & decision$Betweenness_value > betweenness_median
decision$Closeness_above_median <- decision$Valid_metrics & decision$Closeness_value > closeness_median
decision$Degree_above_median <- decision$Valid_metrics & decision$Degree_value > degree_median
decision$Core_include <- decision$Betweenness_above_median &
  decision$Closeness_above_median &
  decision$Degree_above_median
decision$Core_reason <- ifelse(
  !decision$Valid_metrics,
  "excluded: missing node name or topology metric",
  ifelse(
    decision$Core_include,
    "retained: all three metrics strictly above their medians",
    "excluded: at least one metric is not strictly above its median"
  )
)

core_targets <- unique(decision[decision$Core_include, c(
  "Gene_symbol", "Betweenness_value", "Closeness_value", "Degree_value"
)])
core_targets <- core_targets[order(-core_targets$Degree_value, -core_targets$Betweenness_value, core_targets$Gene_symbol), ]

thresholds <- data.frame(
  Metric = c("Betweenness", "Closeness", "Degree", "Valid_nodes", "Core_targets"),
  Threshold_or_count = c(
    betweenness_median,
    closeness_median,
    degree_median,
    sum(decision$Valid_metrics),
    nrow(core_targets)
  ),
  Comparison = c("> median", "> median", "> median", "count", "count"),
  stringsAsFactors = FALSE
)

write.csv(decision, file.path(output_dir, "02_Cytoscape_node_metrics.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(thresholds, file.path(output_dir, "03_core_target_thresholds.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(core_targets, file.path(output_dir, "04_core_targets.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(core_targets$Gene_symbol, file.path(output_dir, "04_core_targets.txt"), useBytes = TRUE)

if (nrow(core_targets) == 0) {
  stop("No core targets satisfy the strict AND rule; inspect the saved metrics and thresholds before a second STRING run.")
}

cat("valid_nodes=", sum(decision$Valid_metrics), "\n", sep = "")
cat("betweenness_median=", signif(betweenness_median, 8), "\n", sep = "")
cat("closeness_median=", signif(closeness_median, 8), "\n", sep = "")
cat("degree_median=", signif(degree_median, 8), "\n", sep = "")
cat("core_targets=", nrow(core_targets), "\n", sep = "")
