args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript build_provisional_stage01_inputs.R <legacy_deduplicated.csv> <output_dir>")

input_file <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(input_file, check.names = FALSE, fileEncoding = "UTF-8")
clean <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z)] <- ""
  z
}
norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
collapse_unique <- function(x) paste(sort(unique(clean(x)[clean(x) != ""])), collapse = "|")

raw[] <- lapply(raw, clean)
raw$pubchem_cid_raw[toupper(raw$pubchem_cid_raw) %in% c("", "N/A", "NA", "NULL", "NONE", "-")] <- ""
raw$cas_raw[toupper(raw$cas_raw) %in% c("", "N/A", "NA", "NULL", "NONE", "-")] <- ""

compound_groups <- split(seq_len(nrow(raw)), paste(raw$herb_pinyin, raw$source, raw$source_ingredient_id, sep = "\r"))
compound_raw <- do.call(rbind, lapply(compound_groups, function(idx) data.frame(
  source = raw$source[idx[[1]]],
  source_herb_id = raw$herb_pinyin[idx[[1]]],
  source_ingredient_id = raw$source_ingredient_id[idx[[1]]],
  herb_normalized_cn = raw$herb_normalized_cn[idx[[1]]],
  herb_pinyin = raw$herb_pinyin[idx[[1]]],
  compound_name_raw = collapse_unique(raw$compound_name_raw[idx]),
  ob = NA_real_,
  dl = NA_real_,
  pubchem_cid_raw = collapse_unique(raw$pubchem_cid_raw[idx]),
  cas_raw = collapse_unique(raw$cas_raw[idx]),
  inchi_key_raw = "",
  evidence_type = "legacy compiled workbook; original screening values unavailable",
  source_url = "",
  access_date = "",
  stringsAsFactors = FALSE
)))
row.names(compound_raw) <- NULL

# Recompute keys from the collapsed values exactly as the Stage 01 script will do.
valid_pubchem <- grepl("^[0-9]+$", compound_raw$pubchem_cid_raw)
valid_cas <- grepl("^[0-9]{2,7}-[0-9]{2}-[0-9]$", compound_raw$cas_raw)
identifier_issue <- (compound_raw$pubchem_cid_raw != "" & !valid_pubchem) |
  (compound_raw$cas_raw != "" & !valid_cas)
identifier_blocking <- (!valid_pubchem & compound_raw$pubchem_cid_raw != "") |
  (compound_raw$pubchem_cid_raw == "" & !valid_cas & compound_raw$cas_raw != "")
compound_raw$source_compound_key <- ifelse(
  valid_pubchem, paste0("PUBCHEM:", compound_raw$pubchem_cid_raw),
  ifelse(valid_cas, paste0("CAS:", compound_raw$cas_raw), paste0("NAME:", norm_key(compound_raw$compound_name_raw)))
)

decision_groups <- split(seq_len(nrow(compound_raw)), compound_raw$source_compound_key)
decisions <- do.call(rbind, lapply(decision_groups, function(idx) {
  key <- compound_raw$source_compound_key[idx[[1]]]
  data.frame(
    source_compound_key = key,
    compound_key = key,
    compound_name_preferred = sort(unique(compound_raw$compound_name_raw[idx]))[[1]],
    name_authority = "legacy workbook; preferred-name authority pending",
    pubchem_cid = collapse_unique(compound_raw$pubchem_cid_raw[idx]),
    cas = collapse_unique(compound_raw$cas_raw[idx]),
    inchi_key = "",
    identity_match_method = ifelse(any(identifier_blocking[idx]), "multiple or malformed identifier pending review",
      ifelse(grepl("^PUBCHEM:", key), "PubChem CID", ifelse(grepl("^CAS:", key), "CAS", "name fallback pending review"))),
    identity_match_status = ifelse(any(identifier_blocking[idx]), "ambiguous",
      ifelse(grepl("^(PUBCHEM|CAS):", key), "exact", "pending")),
    include_compound = TRUE,
    include_reason = "legacy prefiltered import for pipeline test only; raw OB/DL or source evidence must be recovered before freeze",
    stringsAsFactors = FALSE
  )
}))
row.names(decisions) <- NULL

target_raw <- data.frame(
  source = raw$source,
  source_herb_id = raw$herb_pinyin,
  source_ingredient_id = raw$source_ingredient_id,
  compound_key = "",
  target_name_raw = toupper(raw$gene_symbol_input),
  target_id_raw = "",
  evidence_type = "legacy compiled workbook",
  evidence_reference = "",
  herb_selection_tier = "",
  source_url = "",
  access_date = "",
  stringsAsFactors = FALSE
)
target_raw <- unique(target_raw)

compound_raw$source_compound_key <- NULL
write.csv(compound_raw, file.path(output_dir, "compound_raw.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decisions, file.path(output_dir, "compound_decisions_pending_dnp.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(target_raw, file.path(output_dir, "target_raw.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("compound_source_rows=", nrow(compound_raw), "\n", sep = "")
cat("canonical_candidate_keys=", nrow(decisions), "\n", sep = "")
cat("target_source_rows=", nrow(target_raw), "\n", sep = "")
