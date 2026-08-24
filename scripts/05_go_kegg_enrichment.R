if (.Platform$OS.type == "windows" && !grepl("utf8|UTF-8", Sys.getlocale("LC_CTYPE"))) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop(paste(
  "Usage: Rscript 05_go_kegg_enrichment.R",
  "<core_targets.csv|tsv|txt> <background.csv|tsv|txt|NONE> <output_dir>",
  "[compound_target_table|NONE] [theme_keywords.csv|NONE]",
  "[p_cutoff=0.05] [padj_cutoff=0.05] [go_top_n=15] [kegg_top_n=30]",
  "[go_pool_per_ontology=200] [online_timeout_seconds=60]",
  "[input_set_label=target set] [kegg_snapshot_dir|NONE] [kegg_fallback_mode=rest_snapshot|orgdb_path]"
))

required <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

options(stringsAsFactors = FALSE)
core_file <- args[[1]]; background_file <- args[[2]]; output_dir <- args[[3]]
compound_file <- if (length(args) >= 4) args[[4]] else "NONE"
theme_file <- if (length(args) >= 5) args[[5]] else "NONE"
p_cutoff <- if (length(args) >= 6) as.numeric(args[[6]]) else 0.05
padj_cutoff <- if (length(args) >= 7) as.numeric(args[[7]]) else 0.05
go_top_n <- if (length(args) >= 8) as.integer(args[[8]]) else 15L
kegg_top_n <- if (length(args) >= 9) as.integer(args[[9]]) else 30L
go_pool_per_ontology <- if (length(args) >= 10) as.integer(args[[10]]) else 200L
online_timeout_seconds <- if (length(args) >= 11) as.integer(args[[11]]) else 60L
input_set_label <- if (length(args) >= 12) trimws(args[[12]]) else "unspecified target set"
kegg_snapshot_dir <- if (length(args) >= 13) trimws(args[[13]]) else "NONE"
kegg_fallback_mode <- if (length(args) >= 14) trimws(tolower(args[[14]])) else "rest_snapshot"
if (!is.finite(p_cutoff) || p_cutoff <= 0 || p_cutoff > 1) stop("p_cutoff must be in (0, 1].")
if (!is.finite(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) stop("padj_cutoff must be in (0, 1].")
if (go_top_n < 1 || kegg_top_n < 1 || go_pool_per_ontology < go_top_n) stop("Invalid Top-N or GO candidate-pool setting.")
if (online_timeout_seconds < 10) stop("online_timeout_seconds must be at least 10.")
if (!nzchar(input_set_label)) stop("input_set_label must not be empty.")
if (!kegg_fallback_mode %in% c("rest_snapshot", "orgdb_path")) stop("kegg_fallback_mode must be rest_snapshot or orgdb_path.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (ext %in% c("tsv", "txt")) {
    first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
    if (ext == "txt" && length(first) && !grepl("\t", first, fixed = TRUE)) return(data.frame(Gene_symbol = readLines(path, warn = FALSE, encoding = "UTF-8")))
    return(read.delim(path, check.names = FALSE, fileEncoding = "UTF-8"))
  }
  stop("Only CSV, TSV, and TXT inputs are supported: ", path)
}
write_out <- function(x, name) write.csv(x, file.path(output_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
key <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(as.character(x), to = "ASCII//TRANSLIT")))
find_col <- function(df, candidates, required = TRUE) {
  hit <- match(key(candidates), key(names(df)), nomatch = 0); hit <- hit[hit > 0]
  if (!length(hit)) { if (required) stop("Missing column; expected one of: ", paste(candidates, collapse = ", ")); return(NA_character_) }
  names(df)[hit[[1]]]
}
extract_symbols <- function(df) {
  nm <- find_col(df, c("Gene_symbol", "Gene Symbol", "gene_symbol_primary", "symbol", "gene", "preferredName", "name"))
  x <- toupper(trimws(as.character(df[[nm]])))
  sort(unique(x[!is.na(x) & x != "" & x != "NA"]))
}
map_symbols <- function(symbols) {
  out <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
  out <- unique(out[!is.na(out$ENTREZID) & out$ENTREZID != "", c("SYMBOL", "ENTREZID")])
  out[order(out$SYMBOL, out$ENTREZID), , drop = FALSE]
}
append_symbols <- function(df) {
  if (!nrow(df)) { df$Gene_symbols <- character(); return(df) }
  df$geneID <- vapply(strsplit(as.character(df$geneID), "/", fixed = TRUE), function(z) paste(sort(unique(z[z != ""])), collapse = "/"), character(1))
  ids <- unique(unlist(strsplit(as.character(df$geneID), "/", fixed = TRUE)))
  mp <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ids, keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))
  mp <- unique(mp[!is.na(mp$SYMBOL), ]); lookup <- split(mp$SYMBOL, mp$ENTREZID)
  df$Gene_symbols <- vapply(strsplit(as.character(df$geneID), "/", fixed = TRUE), function(z) paste(sort(unique(unlist(lookup[z], use.names = FALSE))), collapse = "/"), character(1))
  df
}
parse_ratio <- function(x) vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) if (length(z) == 2) as.numeric(z[1]) / as.numeric(z[2]) else NA_real_, numeric(1))
rank_terms <- function(df) {
  if (!nrow(df)) return(df)
  df$FDR_significant <- is.finite(df$p.adjust) & df$p.adjust < padj_cutoff
  df$GeneRatio_numeric <- parse_ratio(df$GeneRatio)
  df[order(df$p.adjust, df$pvalue, -df$Count, df$ID, na.last = TRUE), , drop = FALSE]
}
jaccard_reduce <- function(df, cutoff = 0.7) {
  if (!nrow(df)) return(list(kept = df, removed = transform(df, Redundant_with_ID = character(), Jaccard_overlap = numeric(), Removal_reason = character())))
  sets <- strsplit(as.character(df$geneID), "/", fixed = TRUE); kept <- integer(); removed <- list()
  for (i in seq_len(nrow(df))) {
    if (!length(kept)) { kept <- i; next }
    sims <- vapply(kept, function(j) length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]])), numeric(1))
    if (max(sims, na.rm = TRUE) >= cutoff) {
      j <- kept[[which.max(sims)]]; row <- df[i, , drop = FALSE]
      row$Redundant_with_ID <- df$ID[[j]]; row$Jaccard_overlap <- max(sims, na.rm = TRUE); row$Removal_reason <- paste0("gene-set Jaccard overlap >= ", cutoff)
      removed[[length(removed) + 1L]] <- row
    } else kept <- c(kept, i)
  }
  empty <- transform(df[0, , drop = FALSE], Redundant_with_ID = character(), Jaccard_overlap = numeric(), Removal_reason = character())
  list(kept = df[kept, , drop = FALSE], removed = if (length(removed)) do.call(rbind, removed) else empty)
}
balanced_go <- function(df, n) {
  if (!nrow(df)) return(df)
  quota <- floor(n / 3); chosen <- unlist(lapply(c("BP", "MF", "CC"), function(z) head(which(df$ONTOLOGY == z), quota)), use.names = FALSE)
  chosen <- c(chosen, head(setdiff(seq_len(nrow(df)), chosen), n - length(chosen))); df[chosen, , drop = FALSE]
}

core_symbols <- extract_symbols(read_table(core_file)); if (!length(core_symbols)) stop("No valid core-target gene symbols were found.")
core_map <- map_symbols(core_symbols); if (!nrow(core_map)) stop("No core-target symbols mapped to human Entrez IDs.")
write_out(core_map, "01_core_symbol_entrez_mapping.csv")
write_out(data.frame(Gene_symbol = setdiff(core_symbols, core_map$SYMBOL)), "02_core_unmapped_symbols.csv")

background_ids <- NULL; background_mode <- "clusterProfiler default annotated universe (senior workflow: universe = NULL)"
if (toupper(background_file) != "NONE") {
  background_map <- map_symbols(extract_symbols(read_table(background_file)))
  missing_core <- setdiff(unique(core_map$ENTREZID), unique(background_map$ENTREZID))
  if (length(missing_core)) { write_out(core_map[core_map$ENTREZID %in% missing_core, ], "03_core_ids_missing_from_background.csv"); stop("Mapped core targets are absent from the supplied background.") }
  background_ids <- unique(background_map$ENTREZID); background_mode <- "user-supplied background"
  write_out(background_map, "03_background_symbol_entrez_mapping.csv")
} else write_out(data.frame(SYMBOL = character(), ENTREZID = character()), "03_background_symbol_entrez_mapping.csv")

core_ids <- unique(core_map$ENTREZID)
go_obj <- clusterProfiler::enrichGO(gene = core_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", universe = background_ids, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
go_all <- rank_terms(append_symbols(as.data.frame(go_obj)))
write_out(go_all, "04_GO_all_results.csv")
go_sig <- go_all[is.finite(go_all$pvalue) & go_all$pvalue < p_cutoff, , drop = FALSE]; write_out(go_sig, "05_GO_p_lt_cutoff.csv")
go_pool <- do.call(rbind, lapply(c("BP", "MF", "CC"), function(z) head(go_sig[go_sig$ONTOLOGY == z, , drop = FALSE], go_pool_per_ontology)))
outside <- go_sig[!go_sig$ID %in% go_pool$ID, , drop = FALSE]; if (nrow(outside)) outside$Pool_status <- paste0("below top ", go_pool_per_ontology, " redundancy pool in ontology")
go_red <- jaccard_reduce(go_pool); go_kept <- rank_terms(go_red$kept); go_removed <- go_red$removed
write_out(go_removed, "06_GO_redundancy_removed.csv"); write_out(outside, "06b_GO_outside_redundancy_candidate_pool.csv"); write_out(go_kept, "07_GO_nonredundant.csv")
go_top <- balanced_go(go_kept, min(go_top_n, nrow(go_kept))); if (nrow(go_top)) go_top$Display_rank <- seq_len(nrow(go_top))
write_out(head(go_kept, go_top_n), "08_GO_top_overall_ranking.csv"); write_out(go_top, "09_GO_top15_balanced_display.csv")

kegg_error <- ""; kegg_online_error <- ""; kegg_method <- "clusterProfiler::enrichKEGG, organism=hsa"
old_timeout <- getOption("timeout"); options(timeout = online_timeout_seconds)
kegg_obj <- tryCatch(clusterProfiler::enrichKEGG(gene = core_ids, organism = "hsa", universe = background_ids, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"), error = identity)
options(timeout = old_timeout)
if (inherits(kegg_obj, "error")) {
  kegg_online_error <- conditionMessage(kegg_obj)
  kegg_error <- kegg_online_error
  kegg_all <- go_all[0, setdiff(names(go_all), "ONTOLOGY"), drop = FALSE]
} else {
  kegg_raw <- as.data.frame(kegg_obj)
  if (!nrow(kegg_raw)) {
    kegg_online_error <- "KEGG returned no pathway rows; inspect retrieval warnings."
    kegg_error <- kegg_online_error
    kegg_all <- go_all[0, setdiff(names(go_all), "ONTOLOGY"), drop = FALSE]
  } else {
    kegg_all <- rank_terms(append_symbols(kegg_raw))
  }
}
if (!nrow(kegg_all) && toupper(kegg_snapshot_dir) != "NONE") {
  link_file <- file.path(kegg_snapshot_dir, "link_pathway_hsa.tsv")
  name_file <- file.path(kegg_snapshot_dir, "list_pathway_hsa.tsv")
  required_snapshot_files <- if (kegg_fallback_mode == "rest_snapshot") c(link_file, name_file) else c(name_file)
  if (any(!file.exists(required_snapshot_files))) {
    stop("Explicit KEGG fallback is missing required snapshot file(s) in: ", kegg_snapshot_dir)
  }
  names_snapshot <- read.delim(name_file, header = FALSE, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE)
  if (kegg_fallback_mode == "rest_snapshot") {
    links <- read.delim(link_file, header = FALSE, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE)
    term2gene <- unique(data.frame(
      term = sub("^path:", "", as.character(links[[2]])),
      gene = sub("^hsa:", "", as.character(links[[1]])),
      stringsAsFactors = FALSE
    ))
  } else {
    all_entrez_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")
    orgdb_path <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = all_entrez_keys,
      keytype = "ENTREZID",
      columns = "PATH"
    ))
    orgdb_path <- orgdb_path[!is.na(orgdb_path$PATH) & orgdb_path$PATH != "", , drop = FALSE]
    term2gene <- unique(data.frame(
      term = paste0("hsa", as.character(orgdb_path$PATH)),
      gene = as.character(orgdb_path$ENTREZID),
      stringsAsFactors = FALSE
    ))
  }
  term2name <- unique(data.frame(
    term = sub("^path:", "", as.character(names_snapshot[[1]])),
    name = sub(" - Homo sapiens \\(human\\)$", "", as.character(names_snapshot[[2]])),
    stringsAsFactors = FALSE
  ))
  snapshot_obj <- clusterProfiler::enricher(
    gene = core_ids,
    universe = unique(term2gene$gene),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  snapshot_raw <- as.data.frame(snapshot_obj)
  if (nrow(snapshot_raw)) {
    kegg_all <- rank_terms(append_symbols(snapshot_raw))
    kegg_method <- if (kegg_fallback_mode == "rest_snapshot") {
      "clusterProfiler::enricher with frozen official KEGG REST pathway-gene snapshot"
    } else {
      paste0("clusterProfiler::enricher with org.Hs.eg.db PATH mapping v", as.character(utils::packageVersion("org.Hs.eg.db")), " and frozen KEGG pathway names")
    }
    kegg_error <- ""
  } else {
    kegg_error <- "Online enrichKEGG failed and the explicit KEGG snapshot fallback returned no pathways."
  }
}
write_out(kegg_all, "10_KEGG_all_results.csv")
kegg_sig <- kegg_all[is.finite(kegg_all$pvalue) & kegg_all$pvalue < p_cutoff, , drop = FALSE]; write_out(kegg_sig, "11_KEGG_p_lt_cutoff.csv")
kegg_red <- jaccard_reduce(kegg_sig); kegg_kept <- rank_terms(kegg_red$kept); kegg_removed <- kegg_red$removed
kegg_top <- head(kegg_kept, kegg_top_n); if (nrow(kegg_top)) kegg_top$Display_rank <- seq_len(nrow(kegg_top))
write_out(kegg_removed, "12_KEGG_redundancy_removed.csv"); write_out(kegg_kept, "13_KEGG_nonredundant.csv"); write_out(kegg_top, "14_KEGG_top30_display.csv")

wrap_labels <- function(x, width = 48) vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
save_plots <- function(df, prefix, facet = FALSE) {
  if (!nrow(df)) return(FALSE)
  d <- df; d$Description_plot <- factor(wrap_labels(d$Description), levels = rev(wrap_labels(d$Description))); d$minus_log10_p <- -log10(pmax(d$pvalue, .Machine$double.xmin))
  p1 <- ggplot2::ggplot(d, ggplot2::aes(GeneRatio_numeric, Description_plot, size = Count, color = minus_log10_p)) + ggplot2::geom_point() + ggplot2::scale_color_viridis_c() + ggplot2::labs(x = "Gene ratio", y = NULL, color = "-log10(P)") + ggplot2::theme_bw()
  p2 <- ggplot2::ggplot(d, ggplot2::aes(minus_log10_p, Description_plot, fill = GeneRatio_numeric)) + ggplot2::geom_col() + ggplot2::scale_fill_viridis_c() + ggplot2::labs(x = "-log10(P)", y = NULL, fill = "Gene ratio") + ggplot2::theme_bw()
  if (facet) { p1 <- p1 + ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y"); p2 <- p2 + ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") }
  h <- max(6, min(18, 0.32 * nrow(d) + 2)); ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_dot.png")), p1, width = 10, height = h, dpi = 600); ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_dot.pdf")), p1, width = 10, height = h)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_bar.png")), p2, width = 10, height = h, dpi = 600); ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_bar.pdf")), p2, width = 10, height = h); TRUE
}
save_plots(go_top, "15_GO_top15", TRUE); save_plots(kegg_top, "16_KEGG_top30", FALSE)

theme_hits <- data.frame()
if (toupper(theme_file) != "NONE") {
  themes <- read_table(theme_file); theme_col <- find_col(themes, "Theme"); keyword_col <- find_col(themes, "Keyword")
  common <- c("ID", "Description", "pvalue", "p.adjust", "Gene_symbols")
  tag_database <- function(df, label) {
    out <- df[, common, drop = FALSE]
    out$Database <- rep(label, nrow(out))
    out
  }
  combined <- rbind(tag_database(go_kept, "GO"), tag_database(kegg_kept, "KEGG")); hits <- list()
  for (i in seq_len(nrow(themes))) { m <- grepl(as.character(themes[[keyword_col]][i]), combined$Description, ignore.case = TRUE, fixed = TRUE); if (any(m)) { z <- combined[m, ]; z$Theme <- themes[[theme_col]][i]; z$Matched_keyword <- themes[[keyword_col]][i]; hits[[length(hits) + 1L]] <- z } }
  if (length(hits)) theme_hits <- do.call(rbind, hits)
}
write_out(theme_hits, "17_pharmacology_theme_candidates.csv")
compound_core <- data.frame()
if (toupper(compound_file) != "NONE") { compounds <- read_table(compound_file); nm <- find_col(compounds, c("Gene_symbol", "Gene Symbol", "gene_symbol_primary", "symbol", "target")); compounds$Core_Gene_symbol <- toupper(trimws(as.character(compounds[[nm]]))); compound_core <- compounds[compounds$Core_Gene_symbol %in% core_symbols, , drop = FALSE] }
write_out(compound_core, "18_active_compound_core_target_links.csv")

qc <- data.frame(Item = c("Input target set", "Input symbols", "Mapped symbols", "Background mode", "GO all terms", "GO raw-P eligible terms", "GO displayed terms", "KEGG all pathways", "KEGG raw-P eligible pathways", "KEGG displayed pathways", "KEGG online error", "KEGG final error", "GO method", "KEGG method"), Value = c(input_set_label, length(core_symbols), nrow(core_map), background_mode, nrow(go_all), nrow(go_sig), nrow(go_top), nrow(kegg_all), nrow(kegg_sig), nrow(kegg_top), kegg_online_error, kegg_error, "clusterProfiler::enrichGO, ont=ALL", kegg_method))
write_out(qc, "19_enrichment_qc.csv")
parameters <- data.frame(Parameter = c("Input target set", "GO method", "KEGG method", "KEGG fallback mode", "KEGG snapshot directory", "Organism", "Background", "P cutoff", "Adjustment", "GO display", "KEGG display", "Figure resolution"), Value = c(input_set_label, "clusterProfiler::enrichGO, ont=ALL", kegg_method, kegg_fallback_mode, kegg_snapshot_dir, "Homo sapiens", background_mode, p_cutoff, "BH", paste0("balanced Top ", go_top_n), paste0("Top ", kegg_top_n), "PNG 600 dpi and PDF"))
write_out(parameters, "20_enrichment_parameters.csv"); capture.output(sessionInfo(), file = file.path(output_dir, "21_sessionInfo.txt"))
cat("input_target_set=", input_set_label, "\n", sep = ""); cat("input_symbols=", length(core_symbols), "\n", sep = ""); cat("go_displayed=", nrow(go_top), "\n", sep = ""); cat("kegg_displayed=", nrow(kegg_top), "\n", sep = ""); if (nzchar(kegg_error)) cat("kegg_status=failed; see QC table\n")
