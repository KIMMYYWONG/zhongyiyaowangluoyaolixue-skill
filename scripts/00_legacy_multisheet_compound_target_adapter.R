args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 00_legacy_multisheet_compound_target_adapter.R <input.xlsx> <sheet_herb_map.csv> <output.csv>")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) stop("The openxlsx package is required.")

input_file <- args[[1]]
map_file <- args[[2]]
output_file <- args[[3]]
if (!file.exists(input_file)) stop("Workbook not found: ", input_file)
if (!file.exists(map_file)) stop("Sheet-herb map not found: ", map_file)

mapping <- read.csv(map_file, check.names = FALSE, fileEncoding = "UTF-8")
required <- c("sheet_index", "Herb_name", "Herb_code")
if (!all(required %in% names(mapping))) stop("Sheet-herb map must contain: ", paste(required, collapse = ", "))
mapping$sheet_index <- suppressWarnings(as.integer(mapping$sheet_index))
if (any(!is.finite(mapping$sheet_index)) || anyDuplicated(mapping$sheet_index)) stop("sheet_index values must be unique positive integers.")
if (any(is.na(mapping$Herb_name)) || any(trimws(mapping$Herb_name) == "") || anyDuplicated(trimws(mapping$Herb_name))) stop("Herb_name values must be non-empty and unique.")

sheets <- openxlsx::getSheetNames(input_file)
if (any(mapping$sheet_index < 1L | mapping$sheet_index > length(sheets))) stop("Sheet map contains an index outside the workbook sheet range.")

clean <- function(x) {
  y <- trimws(gsub("[[:space:]]+", " ", as.character(x)))
  y[is.na(y)] <- ""
  y
}

records <- list()
for (i in seq_len(nrow(mapping))) {
  sheet_index <- mapping$sheet_index[[i]]
  x <- openxlsx::read.xlsx(input_file, sheet = sheet_index, check.names = FALSE, detectDates = FALSE)
  if (!nrow(x)) next
  if (ncol(x) < 3L) stop("Each legacy sheet must contain at least three columns (ingredient ID, ingredient name, gene): sheet index ", sheet_index)
  z <- data.frame(
    Source_row_id = sprintf("sheet%02d-row%06d", sheet_index, seq_len(nrow(x))),
    Herb_name = clean(mapping$Herb_name[[i]]),
    Herb_code = tolower(clean(mapping$Herb_code[[i]])),
    Ingredient_ID = clean(x[[1]]),
    Ingredient_name = clean(x[[2]]),
    Gene_symbol = toupper(clean(x[[3]])),
    Source = "Legacy imported multi-sheet compound-target workbook",
    Source_sheet_index = sheet_index,
    Source_sheet_name = sheets[[sheet_index]],
    stringsAsFactors = FALSE
  )
  z <- z[z$Ingredient_ID != "" & z$Ingredient_name != "" & z$Gene_symbol != "", , drop = FALSE]
  records[[length(records) + 1L]] <- z
}

if (!length(records)) stop("No valid compound-target rows were extracted.")
out <- do.call(rbind, records)
out <- out[grepl("^[A-Z0-9][A-Z0-9._-]*$", out$Gene_symbol), , drop = FALSE]
out <- out[!duplicated(out[, c("Herb_name", "Ingredient_ID", "Ingredient_name", "Gene_symbol")]), , drop = FALSE]
out <- out[order(out$Herb_name, out$Ingredient_ID, out$Gene_symbol), ]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("sheets=", length(unique(out$Source_sheet_index)), "\n", sep = "")
cat("rows=", nrow(out), "\n", sep = "")
cat("unique_genes=", length(unique(out$Gene_symbol)), "\n", sep = "")
