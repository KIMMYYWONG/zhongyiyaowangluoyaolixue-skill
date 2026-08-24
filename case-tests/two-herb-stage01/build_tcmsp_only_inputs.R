args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript build_tcmsp_only_inputs.R <legacy.csv> <evidence.csv> <uniprot.csv> <output_dir>")

legacy <- read.csv(args[[1]], check.names = FALSE, fileEncoding = "UTF-8")
evidence <- read.csv(args[[2]], check.names = FALSE, fileEncoding = "UTF-8-BOM")
uniprot <- read.csv(args[[3]], check.names = FALSE, fileEncoding = "UTF-8-BOM")
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clean <- function(x) { z <- trimws(as.character(x)); z[is.na(z)] <- ""; z }
legacy[] <- lapply(legacy, clean)
evidence[] <- lapply(evidence, clean)
legacy <- legacy[legacy$source == "TCMSP", ]

key <- function(herb, id) paste(herb, id, sep = "\r")
m <- match(key(legacy$herb_pinyin, legacy$source_ingredient_id), key(evidence$herb_pinyin, evidence$source_ingredient_id))
if (anyNA(m)) stop("TCMSP legacy rows without evidence: ", paste(unique(legacy$source_ingredient_id[is.na(m)]), collapse = ";"))

compound_raw <- unique(data.frame(
  source = "TCMSP",
  source_herb_id = evidence$herb_pinyin,
  source_ingredient_id = evidence$source_ingredient_id,
  herb_normalized_cn = evidence$herb_pinyin,
  herb_pinyin = evidence$herb_pinyin,
  compound_name_raw = evidence$compound_name_preferred,
  ob = as.numeric(evidence$ob),
  dl = as.numeric(evidence$dl),
  pubchem_cid_raw = evidence$pubchem_cid,
  cas_raw = evidence$cas,
  inchi_key_raw = evidence$inchi_key,
  evidence_type = "TCMSP molecule page with recorded cross-check",
  source_url = evidence$tcmsp_url,
  access_date = "2026-08-24",
  stringsAsFactors = FALSE
))

source_key <- ifelse(compound_raw$inchi_key_raw != "", paste0("INCHIKEY:", compound_raw$inchi_key_raw),
  ifelse(compound_raw$pubchem_cid_raw != "", paste0("PUBCHEM:", compound_raw$pubchem_cid_raw),
    paste0("CAS:", compound_raw$cas_raw)))
decision_groups <- split(seq_len(nrow(compound_raw)), source_key)
decisions <- do.call(rbind, lapply(decision_groups, function(idx) data.frame(
  source_compound_key = source_key[idx[[1]]],
  compound_key = source_key[idx[[1]]],
  compound_name_preferred = evidence$compound_name_preferred[idx[[1]]],
  name_authority = "TCMSP/PubChem provisional name; DNP review pending",
  pubchem_cid = paste(sort(unique(evidence$pubchem_cid[idx][evidence$pubchem_cid[idx] != ""])), collapse = "|"),
  cas = paste(sort(unique(evidence$cas[idx][evidence$cas[idx] != ""])), collapse = "|"),
  inchi_key = paste(sort(unique(evidence$inchi_key[idx][evidence$inchi_key[idx] != ""])), collapse = "|"),
  identity_match_method = paste(sort(unique(evidence$identity_match_method[idx])), collapse = "|"),
  identity_match_status = "exact",
  include_compound = TRUE,
  include_reason = "TCMSP OB >= 30 and DL >= 0.18; stable identity reviewed; DNP preferred name pending",
  stringsAsFactors = FALSE
)))
row.names(decisions) <- NULL

target_raw <- unique(data.frame(
  source = "TCMSP",
  source_herb_id = legacy$herb_pinyin,
  source_ingredient_id = legacy$source_ingredient_id,
  compound_key = "",
  target_name_raw = toupper(legacy$gene_symbol_input),
  target_id_raw = "",
  evidence_type = "legacy TCMSP compound-target relationship",
  evidence_reference = "",
  herb_selection_tier = "",
  source_url = "",
  access_date = "",
  stringsAsFactors = FALSE
))

write.csv(compound_raw, file.path(output_dir, "compound_raw.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decisions, file.path(output_dir, "compound_decisions_dnp_pending.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(target_raw, file.path(output_dir, "target_raw.csv"), row.names = FALSE, fileEncoding = "UTF-8")

accepted <- uniprot[tolower(uniprot$mapping_status) == "exact" & uniprot$uniprot_accession != "", ]
tm <- match(target_raw$target_name_raw, accepted$target_name_raw)
provisional <- target_raw[!is.na(tm), ]
evm <- match(key(provisional$source_herb_id, provisional$source_ingredient_id), key(evidence$herb_pinyin, evidence$source_ingredient_id))
um <- match(provisional$target_name_raw, accepted$target_name_raw)
provisional_out <- unique(data.frame(
  herb_pinyin = provisional$source_herb_id,
  tcmsp_molid = provisional$source_ingredient_id,
  compound_name_provisional = evidence$compound_name_preferred[evm],
  compound_name_status = "DNP review pending - not frozen",
  ob = as.numeric(evidence$ob[evm]),
  dl = as.numeric(evidence$dl[evm]),
  pubchem_cid = evidence$pubchem_cid[evm],
  cas = evidence$cas[evm],
  inchi_key = evidence$inchi_key[evm],
  gene_symbol_primary = accepted$gene_symbol_primary[um],
  uniprot_accession = accepted$uniprot_accession[um],
  source = "TCMSP",
  status = "provisional: DNP naming evidence pending",
  source_url = evidence$tcmsp_url[evm],
  stringsAsFactors = FALSE
))
provisional_out <- provisional_out[order(provisional_out$herb_pinyin, provisional_out$tcmsp_molid, provisional_out$gene_symbol_primary), ]
write.csv(provisional_out, file.path(output_dir, "provisional_tcmsp_compound_target.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summary <- do.call(rbind, lapply(split(provisional_out, provisional_out$herb_pinyin), function(x) data.frame(
  herb_pinyin = x$herb_pinyin[[1]],
  retained_tcmsp_compounds = length(unique(x$tcmsp_molid)),
  provisional_compound_target_rows = nrow(x),
  unique_uniprot_targets = length(unique(x$gene_symbol_primary)),
  stringsAsFactors = FALSE
)))
write.csv(summary, file.path(output_dir, "provisional_tcmsp_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("tcmsp_compounds=", nrow(compound_raw), "\n", sep = "")
cat("provisional_rows=", nrow(provisional_out), "\n", sep = "")
