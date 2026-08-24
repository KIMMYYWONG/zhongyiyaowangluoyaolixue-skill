if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop(paste(
  "Usage: Rscript 01_drug_target_normalize_filter.R",
  "<compound_raw.csv|tsv> <target_raw.csv|tsv> <compound_decisions.csv|tsv|NONE>",
  "<uniprot_mapping.csv|tsv> <output_dir> [ob_min=30] [dl_min=0.18]"
))
options(stringsAsFactors = FALSE)
compound_file <- args[[1]]; target_file <- args[[2]]; decision_file <- args[[3]]
mapping_file <- args[[4]]; output_dir <- args[[5]]
ob_min <- if (length(args) >= 6) as.numeric(args[[6]]) else 30
dl_min <- if (length(args) >= 7) as.numeric(args[[7]]) else 0.18
if (!is.finite(ob_min) || !is.finite(dl_min)) stop("OB and DL thresholds must be numeric.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (ext %in% c("tsv", "txt")) return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  stop("Only CSV, TSV, and TXT inputs are supported: ", path)
}
write_out <- function(x, name) write.csv(x, file.path(output_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
find_col <- function(df, candidates, required = TRUE) {
  hit <- match(norm_key(candidates), norm_key(names(df)), nomatch = 0); hit <- hit[hit > 0]
  if (!length(hit)) { if (required) stop("Missing column; expected one of: ", paste(candidates, collapse = ", ")); return(NA_character_) }
  names(df)[hit[[1]]]
}
clean <- function(x) trimws(gsub("[[:space:]]+", " ", as.character(x)))
blank <- function(n) rep("", n)
get_col <- function(df, candidates, default = "", required = FALSE) {
  nm <- find_col(df, candidates, required = required)
  if (is.na(nm)) return(rep(default, nrow(df)))
  value <- clean(df[[nm]])
  value[is.na(value)] <- default
  value
}
as_bool <- function(x) tolower(clean(x)) %in% c("true", "t", "1", "yes", "y", "include", "included")
source_class <- function(x) {
  z <- toupper(gsub("[^A-Z0-9]", "", iconv(clean(x), to = "ASCII//TRANSLIT")))
  ifelse(grepl("TCMSP", z), "TCMSP", ifelse(grepl("SYMMAP", z), "SYMMAP", ifelse(grepl("HERB", z), "HERB2", z)))
}

compound_raw <- read_table(compound_file)
compounds <- data.frame(
  Source_row = seq_len(nrow(compound_raw)),
  source = source_class(get_col(compound_raw, c("source"), required = TRUE)),
  source_herb_id = get_col(compound_raw, c("source_herb_id", "Herb_ID")),
  source_ingredient_id = get_col(compound_raw, c("source_ingredient_id", "Mol_ID", "Ingredient_ID"), required = TRUE),
  herb_normalized_cn = get_col(compound_raw, c("herb_normalized_cn", "Herb_name", "Herb"), required = TRUE),
  herb_pinyin = get_col(compound_raw, c("herb_pinyin", "Herb_pinyin", "Pinyin")),
  compound_name_raw = get_col(compound_raw, c("compound_name_raw", "Compound_name", "Ingredient_name"), required = TRUE),
  ob = suppressWarnings(as.numeric(get_col(compound_raw, c("ob", "OB"), default = NA_character_))),
  dl = suppressWarnings(as.numeric(get_col(compound_raw, c("dl", "DL"), default = NA_character_))),
  pubchem_cid_raw = get_col(compound_raw, c("pubchem_cid_raw", "PubChem_CID", "CID")),
  cas_raw = get_col(compound_raw, c("cas_raw", "CAS", "CAS_number")),
  inchi_key_raw = toupper(get_col(compound_raw, c("inchi_key_raw", "InChIKey"))),
  evidence_type = get_col(compound_raw, c("evidence_type", "Evidence_type")),
  source_url = get_col(compound_raw, c("source_url", "URL")),
  access_date = get_col(compound_raw, c("access_date", "Access_date")),
  stringsAsFactors = FALSE
)
identifier_placeholder <- function(x) toupper(clean(x)) %in% c("", "N/A", "NA", "NULL", "NONE", "-")
compounds$pubchem_cid_raw[identifier_placeholder(compounds$pubchem_cid_raw)] <- ""
compounds$cas_raw[identifier_placeholder(compounds$cas_raw)] <- ""
compounds$inchi_key_raw[identifier_placeholder(compounds$inchi_key_raw)] <- ""
compound_dedupe_cols <- setdiff(names(compounds), "Source_row")
compound_duplicate_flag <- duplicated(compounds[, compound_dedupe_cols, drop = FALSE])
write_out(compounds[compound_duplicate_flag, , drop = FALSE], "00a_compound_exact_duplicates_removed.csv")
compounds <- compounds[!compound_duplicate_flag, , drop = FALSE]
valid_inchi_key <- grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", compounds$inchi_key_raw)
valid_pubchem_cid <- grepl("^[0-9]+$", compounds$pubchem_cid_raw)
valid_cas <- grepl("^[0-9]{2,7}-[0-9]{2}-[0-9]$", compounds$cas_raw)
compounds$source_identifier_issue <-
  (compounds$inchi_key_raw != "" & !valid_inchi_key) |
  (compounds$pubchem_cid_raw != "" & !valid_pubchem_cid) |
  (compounds$cas_raw != "" & !valid_cas)
compounds$source_identifier_blocking <-
  (!valid_inchi_key & compounds$inchi_key_raw != "") |
  (!valid_inchi_key & compounds$inchi_key_raw == "" & !valid_pubchem_cid & compounds$pubchem_cid_raw != "") |
  (!valid_inchi_key & compounds$inchi_key_raw == "" & compounds$pubchem_cid_raw == "" & !valid_cas & compounds$cas_raw != "")
compounds$source_compound_key <- ifelse(valid_inchi_key, paste0("INCHIKEY:", compounds$inchi_key_raw),
  ifelse(valid_pubchem_cid, paste0("PUBCHEM:", compounds$pubchem_cid_raw),
    ifelse(valid_cas, paste0("CAS:", compounds$cas_raw), paste0("NAME:", norm_key(compounds$compound_name_raw)))))
compounds$compound_key <- compounds$source_compound_key
compounds$identity_match_method <- ifelse(grepl("^INCHIKEY:", compounds$compound_key), "InChIKey",
  ifelse(grepl("^PUBCHEM:", compounds$compound_key), "PubChem CID",
    ifelse(grepl("^CAS:", compounds$compound_key), "CAS",
      ifelse(compounds$source_identifier_issue, "multiple or malformed identifier pending review", "name fallback pending review"))))
compounds$identity_match_status <- ifelse(compounds$source_identifier_blocking, "ambiguous",
  ifelse(grepl("^(INCHIKEY|PUBCHEM|CAS):", compounds$compound_key), "exact", "pending"))
valid_compound <- compounds$source != "" & compounds$source_ingredient_id != "" & compounds$herb_normalized_cn != "" & compounds$compound_name_raw != "" & compounds$source_compound_key != "NAME:"
compounds$source_rule <- ifelse(compounds$source == "TCMSP", paste0("OB >= ", ob_min, " AND DL >= ", dl_min),
  ifelse(compounds$source == "SYMMAP", paste0("OB >= ", ob_min), ifelse(compounds$source == "HERB2", "no OB/DL; relationship requires literature evidence", "unsupported source")))
compounds$source_rule_pass <- valid_compound & ifelse(compounds$source == "TCMSP", !is.na(compounds$ob) & compounds$ob >= ob_min & !is.na(compounds$dl) & compounds$dl >= dl_min,
  ifelse(compounds$source == "SYMMAP", !is.na(compounds$ob) & compounds$ob >= ob_min, compounds$source == "HERB2"))
compounds$compound_name_preferred <- compounds$compound_name_raw
compounds$name_authority <- "raw source name pending normalization"
compounds$manual_include <- NA
compounds$manual_reason <- ""

if (toupper(decision_file) != "NONE") {
  decisions <- read_table(decision_file)
  decision_key_col <- find_col(decisions, c("compound_key"))
  source_key_col <- find_col(decisions, c("source_compound_key", "raw_compound_key"), required = FALSE)
  preferred_col <- find_col(decisions, c("compound_name_preferred", "Preferred_name"))
  authority_col <- find_col(decisions, c("name_authority", "Naming_authority"), required = FALSE)
  accepted_pubchem_col <- find_col(decisions, c("pubchem_cid", "accepted_pubchem_cid"), required = FALSE)
  accepted_cas_col <- find_col(decisions, c("cas", "accepted_cas"), required = FALSE)
  accepted_inchi_col <- find_col(decisions, c("inchi_key", "accepted_inchi_key"), required = FALSE)
  identity_method_col <- find_col(decisions, c("identity_match_method", "match_method"), required = FALSE)
  identity_status_col <- find_col(decisions, c("identity_match_status", "match_status"), required = FALSE)
  include_col <- find_col(decisions, c("include_compound", "Include"), required = FALSE)
  reason_col <- find_col(decisions, c("include_reason", "Reason"), required = FALSE)
  d <- data.frame(source_compound_key = if (is.na(source_key_col)) clean(decisions[[decision_key_col]]) else clean(decisions[[source_key_col]]),
    compound_key = clean(decisions[[decision_key_col]]),
    preferred = clean(decisions[[preferred_col]]),
    authority = if (is.na(authority_col)) blank(nrow(decisions)) else clean(decisions[[authority_col]]),
    accepted_pubchem = if (is.na(accepted_pubchem_col)) blank(nrow(decisions)) else clean(decisions[[accepted_pubchem_col]]),
    accepted_cas = if (is.na(accepted_cas_col)) blank(nrow(decisions)) else clean(decisions[[accepted_cas_col]]),
    accepted_inchi = if (is.na(accepted_inchi_col)) blank(nrow(decisions)) else toupper(clean(decisions[[accepted_inchi_col]])),
    identity_method = if (is.na(identity_method_col)) blank(nrow(decisions)) else clean(decisions[[identity_method_col]]),
    identity_status = if (is.na(identity_status_col)) blank(nrow(decisions)) else tolower(clean(decisions[[identity_status_col]])),
    manual_include = if (is.na(include_col)) rep(NA, nrow(decisions)) else as_bool(decisions[[include_col]]),
    manual_reason = if (is.na(reason_col)) blank(nrow(decisions)) else clean(decisions[[reason_col]]))
  if (anyDuplicated(d$source_compound_key)) stop("compound_decisions contains duplicate source_compound_key rows.")
  if (any(d$source_compound_key == "" | d$compound_key == "")) stop("compound_decisions contains an empty source or canonical compound key.")
  m <- match(compounds$source_compound_key, d$source_compound_key)
  has <- !is.na(m)
  matched_rows <- which(has)
  matched_decisions <- m[matched_rows]
  compounds$compound_key[matched_rows] <- d$compound_key[matched_decisions]
  assign_nonblank <- function(destination, values) {
    use <- matched_rows[values != "" & !is.na(values)]
    destination[use] <- values[values != "" & !is.na(values)]
    destination
  }
  compounds$compound_name_preferred <- assign_nonblank(compounds$compound_name_preferred, d$preferred[matched_decisions])
  compounds$name_authority <- assign_nonblank(compounds$name_authority, d$authority[matched_decisions])
  compounds$pubchem_cid_raw <- assign_nonblank(compounds$pubchem_cid_raw, d$accepted_pubchem[matched_decisions])
  compounds$cas_raw <- assign_nonblank(compounds$cas_raw, d$accepted_cas[matched_decisions])
  compounds$inchi_key_raw <- assign_nonblank(compounds$inchi_key_raw, d$accepted_inchi[matched_decisions])
  compounds$identity_match_method <- assign_nonblank(compounds$identity_match_method, d$identity_method[matched_decisions])
  compounds$identity_match_status <- assign_nonblank(compounds$identity_match_status, d$identity_status[matched_decisions])
  compounds$manual_include[matched_rows] <- d$manual_include[matched_decisions]
  compounds$manual_reason[matched_rows] <- d$manual_reason[matched_decisions]
}

collapse_unique <- function(x) {
  z <- sort(unique(clean(x)))
  z <- z[z != "" & !is.na(z)]
  paste(z, collapse = "|")
}
compound_groups <- split(seq_len(nrow(compounds)), compounds$compound_key)
compound_dictionary <- do.call(rbind, lapply(names(compound_groups), function(k) {
  idx <- compound_groups[[k]]
  inchi <- unique(compounds$inchi_key_raw[idx][compounds$inchi_key_raw[idx] != ""])
  pubchem <- unique(compounds$pubchem_cid_raw[idx][compounds$pubchem_cid_raw[idx] != ""])
  cas <- unique(compounds$cas_raw[idx][compounds$cas_raw[idx] != ""])
  preferred <- unique(compounds$compound_name_preferred[idx][compounds$compound_name_preferred[idx] != ""])
  authority <- unique(compounds$name_authority[idx][compounds$name_authority[idx] != ""])
  identity_conflict <- length(inchi) > 1L || length(pubchem) > 1L || length(preferred) > 1L || length(authority) > 1L
  data.frame(compound_key = k, compound_name_preferred = collapse_unique(preferred),
    name_authority = collapse_unique(authority), inchi_key = collapse_unique(inchi),
    pubchem_cid = collapse_unique(pubchem), cas = collapse_unique(cas),
    identity_match_method = collapse_unique(compounds$identity_match_method[idx]),
    identity_match_status = ifelse(identity_conflict, "ambiguous", collapse_unique(compounds$identity_match_status[idx])),
    identity_conflict = identity_conflict,
    source_count = length(unique(compounds$source[idx])), sources = collapse_unique(compounds$source[idx]),
    source_record_count = length(idx), source_ingredient_ids = collapse_unique(compounds$source_ingredient_id[idx]),
    tcmsp_molids = collapse_unique(compounds$source_ingredient_id[idx][compounds$source[idx] == "TCMSP"]),
    herbs = collapse_unique(compounds$herb_normalized_cn[idx]), raw_names = collapse_unique(compounds$compound_name_raw[idx]),
    stringsAsFactors = FALSE)
}))
row.names(compound_dictionary) <- NULL
write_out(compound_dictionary, "01a_compound_identity_dictionary.csv")
write_out(compound_dictionary[compound_dictionary$source_count > 1L | compound_dictionary$identity_conflict, , drop = FALSE],
  "01b_compound_cross_source_dedup_audit.csv")
dict_match <- match(compounds$compound_key, compound_dictionary$compound_key)
compounds$identity_conflict <- compound_dictionary$identity_conflict[dict_match]
compounds$tcmsp_molids <- compound_dictionary$tcmsp_molids[dict_match]
compounds$identity_ready <- grepl("^(INCHIKEY|PUBCHEM|CAS):", compounds$compound_key) &
  compounds$identity_match_status %in% c("exact", "accepted", "reviewed") & !compounds$identity_conflict
compounds$naming_ready <- compounds$compound_name_preferred != "" &
  grepl("dictionary of natural products|天然产物词典|(^|[^a-z])dnp([^a-z]|$)", compounds$name_authority, ignore.case = TRUE) &
  !grepl("pending|not confirmed|unconfirmed|待核验|待确认|未核验|未确认", compounds$name_authority, ignore.case = TRUE)
compounds$screening_candidate <- ifelse(is.na(compounds$manual_include), compounds$source_rule_pass, compounds$manual_include)
compounds$include_compound <- compounds$screening_candidate & compounds$identity_ready & compounds$naming_ready
compounds$include_reason <- ifelse(!compounds$screening_candidate, "excluded by source-specific screening or manual decision",
  ifelse(!compounds$identity_ready, "excluded: unique compound identity unresolved or conflicting",
    ifelse(!compounds$naming_ready, "excluded: Dictionary of Natural Products preferred name not confirmed",
      ifelse(!is.na(compounds$manual_include) & compounds$manual_reason != "", compounds$manual_reason, "retained after screening, identity resolution, and preferred-name review"))))
compounds$manual_override <- !is.na(compounds$manual_include) & compounds$manual_include != compounds$source_rule_pass
write_out(compounds, "01_compound_screening_decisions.csv")

target_raw <- read_table(target_file)
targets <- data.frame(
  Target_source_row = seq_len(nrow(target_raw)),
  source = source_class(get_col(target_raw, c("source"), required = TRUE)),
  source_herb_id = get_col(target_raw, c("source_herb_id", "Herb_ID")),
  source_ingredient_id = get_col(target_raw, c("source_ingredient_id", "Mol_ID", "Ingredient_ID"), required = TRUE),
  compound_key_input = get_col(target_raw, c("compound_key")),
  target_name_raw = get_col(target_raw, c("target_name_raw", "Target_name", "Target"), required = TRUE),
  target_id_raw = get_col(target_raw, c("target_id_raw", "Target_ID")),
  evidence_type = get_col(target_raw, c("evidence_type", "Evidence_type")),
  evidence_reference = get_col(target_raw, c("evidence_reference", "PMID", "DOI", "Reference")),
  herb_selection_tier = get_col(target_raw, c("herb_selection_tier", "HERB_selection_tier")),
  source_url_target = get_col(target_raw, c("source_url", "URL")),
  access_date_target = get_col(target_raw, c("access_date", "Access_date")), stringsAsFactors = FALSE
)
target_dedupe_cols <- setdiff(names(targets), "Target_source_row")
target_duplicate_flag <- duplicated(targets[, target_dedupe_cols, drop = FALSE])
write_out(targets[target_duplicate_flag, , drop = FALSE], "00b_target_exact_duplicates_removed.csv")
targets <- targets[!target_duplicate_flag, , drop = FALSE]
compound_join <- compounds[, c("source", "source_herb_id", "source_ingredient_id", "source_compound_key", "compound_key", "herb_normalized_cn", "herb_pinyin",
  "compound_name_preferred", "name_authority", "pubchem_cid_raw", "cas_raw", "inchi_key_raw", "ob", "dl", "source_rule", "source_rule_pass",
  "identity_match_method", "identity_match_status", "identity_conflict", "tcmsp_molids", "screening_candidate", "include_compound", "include_reason", "manual_override", "source_url", "access_date")]
compound_join <- compound_join[!duplicated(compound_join[, c("source", "source_herb_id", "source_ingredient_id")]), ]
targets <- merge(targets, compound_join, by = c("source", "source_herb_id", "source_ingredient_id"), all.x = TRUE)
targets$compound_key_conflict <- targets$compound_key_input != "" & !is.na(targets$compound_key) &
  targets$compound_key != "" & targets$compound_key_input != targets$source_compound_key & targets$compound_key_input != targets$compound_key
targets$compound_key <- ifelse(!is.na(targets$compound_key) & targets$compound_key != "", targets$compound_key, targets$compound_key_input)

mapping_raw <- read_table(mapping_file)
mapping <- data.frame(
  target_name_raw = clean(mapping_raw[[find_col(mapping_raw, c("target_name_raw", "Target_name", "Target"))]]),
  uniprot_accession = get_col(mapping_raw, c("uniprot_accession", "Entry", "UniProt")),
  gene_symbol_primary = toupper(get_col(mapping_raw, c("gene_symbol_primary", "Gene_symbol", "Gene Names primary"))),
  reviewed_status = get_col(mapping_raw, c("reviewed_status", "Reviewed")),
  organism_result = get_col(mapping_raw, c("organism_result", "Organism")),
  mapping_status = tolower(get_col(mapping_raw, c("mapping_status", "Status"))),
  mapping_note = get_col(mapping_raw, c("mapping_note", "Note")), stringsAsFactors = FALSE
)
accepted_status <- mapping$mapping_status %in% c("exact", "alias", "accepted")
human_result <- grepl("homo sapiens|human|9606", mapping$organism_result, ignore.case = TRUE)
mapping$mapping_accepted <- accepted_status & human_result & mapping$gene_symbol_primary != "" & mapping$uniprot_accession != ""
accepted_multi <- aggregate(gene_symbol_primary ~ target_name_raw, mapping[mapping$mapping_accepted, ], function(x) length(unique(x)))
ambiguous_names <- accepted_multi$target_name_raw[accepted_multi$gene_symbol_primary > 1]
mapping$mapping_accepted[mapping$target_name_raw %in% ambiguous_names] <- FALSE
mapping$mapping_status[mapping$target_name_raw %in% ambiguous_names] <- "ambiguous"
accepted_map <- mapping[mapping$mapping_accepted, ]
accepted_map <- accepted_map[!duplicated(accepted_map$target_name_raw), ]
targets <- merge(targets, accepted_map[, c("target_name_raw", "uniprot_accession", "gene_symbol_primary", "reviewed_status", "organism_result", "mapping_status", "mapping_note")],
  by = "target_name_raw", all.x = TRUE)
targets$mapping_accepted <- !is.na(targets$gene_symbol_primary) & targets$gene_symbol_primary != ""
targets$herb_literature_supported <- targets$source != "HERB2" | targets$evidence_reference != ""
valid_herb_tiers <- c("", "tcmsp_target_overlap", "representative_ingredient_target")
targets$herb_volume_pass <- targets$source != "HERB2" | targets$herb_selection_tier %in% valid_herb_tiers
targets$include_relationship <- !is.na(targets$include_compound) & targets$include_compound & !targets$compound_key_conflict & targets$mapping_accepted &
  targets$herb_literature_supported & targets$herb_volume_pass
targets$relationship_reason <- ifelse(is.na(targets$include_compound), "excluded: compound record not matched",
  ifelse(!targets$include_compound, "excluded: compound not retained",
    ifelse(targets$compound_key_conflict, "excluded: target compound key conflicts with resolved compound identity",
      ifelse(!targets$mapping_accepted, "excluded: target mapping not accepted human mapping",
      ifelse(!targets$herb_literature_supported, "excluded: HERB relationship lacks literature reference",
        ifelse(!targets$herb_volume_pass, "excluded: HERB volume-control tier", "retained"))))))
write_out(targets, "02_target_mapping_and_relationship_decisions.csv")
write_out(mapping[!mapping$mapping_accepted, ], "03_target_mapping_exceptions.csv")

final <- targets[targets$include_relationship, ]
final_out <- unique(data.frame(
  herb_pinyin = final$herb_pinyin,
  pubchem_cid = final$pubchem_cid_raw,
  tcmsp_molid = final$tcmsp_molids,
  compound_name_preferred = final$compound_name_preferred,
  gene_symbol_primary = final$gene_symbol_primary,
  cas = final$cas_raw,
  source = final$source,
  herb_normalized_cn = final$herb_normalized_cn,
  source_herb_id = final$source_herb_id,
  source_ingredient_id = final$source_ingredient_id,
  compound_key = final$compound_key,
  uniprot_accession = final$uniprot_accession,
  evidence_type = final$evidence_type,
  evidence_reference = final$evidence_reference,
  source_url = ifelse(final$source_url_target != "", final$source_url_target, final$source_url),
  access_date = ifelse(final$access_date_target != "", final$access_date_target, final$access_date),
  compound_decision = final$include_reason,
  target_mapping_status = final$mapping_status,
  herb_selection_tier = final$herb_selection_tier,
  stringsAsFactors = FALSE
))
final_out <- final_out[order(final_out$herb_normalized_cn, final_out$compound_name_preferred, final_out$gene_symbol_primary, final_out$source), ]
write_out(final_out, "04_compound_target_final.csv")
dedup_groups <- split(seq_len(nrow(final_out)), paste(final_out$herb_normalized_cn, final_out$compound_key, final_out$gene_symbol_primary, sep = "\r"))
final_dedup <- if (!length(dedup_groups)) data.frame() else do.call(rbind, lapply(dedup_groups, function(idx) data.frame(
  herb_pinyin = collapse_unique(final_out$herb_pinyin[idx]), herb_normalized_cn = collapse_unique(final_out$herb_normalized_cn[idx]),
  compound_key = final_out$compound_key[idx[[1]]], compound_name_preferred = collapse_unique(final_out$compound_name_preferred[idx]),
  gene_symbol_primary = final_out$gene_symbol_primary[idx[[1]]], uniprot_accession = collapse_unique(final_out$uniprot_accession[idx]),
  pubchem_cid = collapse_unique(final_out$pubchem_cid[idx]), tcmsp_molid = collapse_unique(final_out$tcmsp_molid[idx]),
  cas = collapse_unique(final_out$cas[idx]), sources = collapse_unique(final_out$source[idx]),
  source_count = length(unique(final_out$source[idx])), source_ingredient_ids = collapse_unique(final_out$source_ingredient_id[idx]),
  evidence_references = collapse_unique(final_out$evidence_reference[idx]), source_urls = collapse_unique(final_out$source_url[idx]),
  stringsAsFactors = FALSE)))
if (nrow(final_dedup)) row.names(final_dedup) <- NULL
write_out(final_dedup, "04a_compound_target_deduplicated.csv")
drug_targets <- unique(final_out[, c("gene_symbol_primary", "uniprot_accession")])
names(drug_targets)[[1]] <- "Gene_symbol"
drug_targets <- drug_targets[order(drug_targets$Gene_symbol), ]
write_out(drug_targets, "05_unique_drug_targets.csv")

count_one <- function(df, group_cols, label) {
  empty <- data.frame(source = character(), herb_normalized_cn = character(), Count = integer(), Measure = character())
  if (!nrow(df)) return(empty)
  z <- aggregate(rep(1, nrow(df)), df[group_cols], length)
  if (!nrow(z)) return(empty)
  names(z)[ncol(z)] <- "Count"; z$Measure <- label; z
}
counts <- rbind(
  count_one(compounds, c("source", "herb_normalized_cn"), "raw compound rows"),
  count_one(compounds[compounds$include_compound, ], c("source", "herb_normalized_cn"), "retained compound rows"),
  count_one(targets, c("source", "herb_normalized_cn"), "raw target rows"),
  count_one(targets[targets$include_relationship, ], c("source", "herb_normalized_cn"), "retained compound-target rows")
)
write_out(counts, "06_stage01_counts.csv")
manual_evidence_class <- grepl("literature|pharmacopoe|manual addition", compounds$manual_reason, ignore.case = TRUE)
tcmsp_screening_evidence_ok <- compounds$source != "TCMSP" | !compounds$screening_candidate |
  compounds$source_rule_pass | (!is.na(compounds$manual_include) & compounds$manual_include & manual_evidence_class)
symmap_screening_evidence_ok <- compounds$source != "SYMMAP" | !compounds$screening_candidate |
  compounds$source_rule_pass | (!is.na(compounds$manual_include) & compounds$manual_include & manual_evidence_class)
herb_candidate_reference_ok <- targets$source != "HERB2" | is.na(targets$screening_candidate) |
  !targets$screening_candidate | targets$evidence_reference != ""
qc <- data.frame(Check = c("Supported sources only", "TCMSP screening candidates have OB/DL or explicit supplementary evidence",
  "SYMMAP screening candidates have OB or explicit supplementary evidence", "Screening candidates have stable unique compound identifiers",
  "No compound identity conflicts", "Screening candidates have Dictionary of Natural Products preferred names",
  "No target-to-compound key conflicts", "No retained rows with missing gene symbol", "No ambiguous accepted mappings",
  "HERB screening candidates have literature references", "HERB retained rows have literature references",
  "TCMSP rule applied", "SYMMAP rule applied", "Final table nonempty", "Deduplicated relation keys are unique"),
  Pass = c(all(compounds$source %in% c("TCMSP", "SYMMAP", "HERB2")),
    all(tcmsp_screening_evidence_ok), all(symmap_screening_evidence_ok),
    all(!compounds$screening_candidate | compounds$identity_ready), !any(compound_dictionary$identity_conflict),
    all(!compounds$screening_candidate | compounds$naming_ready), !any(targets$compound_key_conflict),
    !any(final_out$gene_symbol_primary == "" | is.na(final_out$gene_symbol_primary)),
    !any(mapping$mapping_accepted & mapping$mapping_status == "ambiguous"),
    all(herb_candidate_reference_ok),
    !any(final_out$source == "HERB2" & final_out$evidence_reference == ""),
    all(compounds$source != "TCMSP" | compounds$source_rule == paste0("OB >= ", ob_min, " AND DL >= ", dl_min)),
    all(compounds$source != "SYMMAP" | compounds$source_rule == paste0("OB >= ", ob_min)), nrow(final_out) > 0,
    !nrow(final_dedup) || !anyDuplicated(final_dedup[, c("herb_normalized_cn", "compound_key", "gene_symbol_primary")])),
  Detail = c("", paste(compounds$source_ingredient_id[!tcmsp_screening_evidence_ok], collapse = ";"),
    paste(compounds$source_ingredient_id[!symmap_screening_evidence_ok], collapse = ";"),
    paste(compounds$source_compound_key[compounds$screening_candidate & !compounds$identity_ready], collapse = ";"),
    paste(compound_dictionary$compound_key[compound_dictionary$identity_conflict], collapse = ";"),
    paste(compounds$compound_key[compounds$screening_candidate & !compounds$naming_ready], collapse = ";"),
    paste(targets$Target_source_row[targets$compound_key_conflict], collapse = ";"), "", paste(ambiguous_names, collapse = ";"),
    paste(targets$Target_source_row[!herb_candidate_reference_ok], collapse = ";"), "",
    paste0("OB=", ob_min, ";DL=", dl_min), paste0("OB=", ob_min), paste0("rows=", nrow(final_out)), paste0("rows=", nrow(final_dedup))))
write_out(qc, "07_stage01_qc.csv")
if (!all(qc$Pass)) stop("Stage 01 QC failed; inspect 07_stage01_qc.csv.")
cat("retained_compounds=", sum(compounds$include_compound), "\n", sep = "")
cat("final_compound_target_rows=", nrow(final_out), "\n", sep = "")
cat("unique_drug_targets=", nrow(drug_targets), "\n", sep = "")
