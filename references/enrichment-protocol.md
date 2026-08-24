# clusterProfiler GO and KEGG enrichment protocol

## Confirmed method

Run both GO and KEGG with `clusterProfiler`, following the confirmed senior workflow. Do not use Metascape unless the user explicitly changes the method.

Run two complete and independent enrichment jobs:

1. the frozen drug-disease intersection target list;
2. the frozen Cytoscape core-target list.

For each job, convert human Gene Symbols to Entrez IDs with `bitr`/`org.Hs.eg.db`, save successful mappings and unmapped symbols, then use that job's unique Entrez IDs for GO and KEGG. Never pool the two input lists, reuse one job's ID-conversion table for the other, or overwrite one job's outputs.

Use separate output directories such as `05A_intersection_GO_KEGG` and `05B_core_GO_KEGG`. Every QC and parameter table must contain the corresponding input-set label.

## Analysis settings

- GO: `clusterProfiler::enrichGO`, `OrgDb = org.Hs.eg.db`, `keyType = ENTREZID`, `ont = ALL`.
- KEGG: `clusterProfiler::enrichKEGG`, `organism = hsa`.
- Multiple-testing adjustment: Benjamini–Hochberg.
- Capture settings: `pvalueCutoff = 1` and `qvalueCutoff = 1`, followed by explicit filtering in saved tables.
- Default background: `universe = NULL`, matching the senior script. A custom background is allowed only when explicitly configured and must contain all mapped input targets.
- Requested display eligibility: raw `P < 0.05`; always retain and report adjusted P and q values.

KEGG may require network access. A retrieval failure must be recorded in QC and must not trigger a silent switch to Metascape or another enrichment source.

When `enrichKEGG` fails because the KEGG REST service is unavailable, an explicitly supplied frozen annotation source may be used for a reproducible retry with `clusterProfiler::enricher`. Two labeled modes are supported: `rest_snapshot` uses frozen official KEGG REST human pathway-gene links and names; `orgdb_path` follows the senior workflow by using the versioned `org.Hs.eg.db` PATH mapping with frozen KEGG pathway names. Preserve retrieval dates, package versions, and hashes, and label the output method rather than calling it `enrichKEGG`. This is an explicit technical fallback within clusterProfiler, not an automatic switch to Metascape or another analysis platform. A run without an explicitly supplied snapshot records the online failure and leaves KEGG outputs empty.

## Redundancy and display

Save complete, raw-P eligible, nonredundant, removed, and display tables separately. Apply display de-duplication after statistical filtering; it does not change the complete clusterProfiler result.

For GO, rank BP, MF, and CC terms by adjusted P, raw P, descending count, and stable ID. Place the best 200 terms per ontology into the default display-candidate pool and apply greedy member-gene Jaccard reduction at `0.7`. Display up to 15 terms, taking up to five each from BP, MF, and CC and filling unused positions by overall rank.

For KEGG, rank by adjusted P, raw P, descending count, and stable ID, then apply greedy member-gene Jaccard reduction at `0.7`. Display up to 30 pathways. If fewer terms qualify, plot all and report the shortfall.

Export dot and bar figures as 600-dpi PNG and vector PDF. Use [pharmacology-theme-keywords.template.csv](../assets/pharmacology-theme-keywords.template.csv) only to retrieve actually enriched terms for interpretation; keywords never force inclusion.

## Runnable module

Run [scripts/05_go_kegg_enrichment.R](../scripts/05_go_kegg_enrichment.R):

```text
Rscript 05_go_kegg_enrichment.R <targets.csv|tsv|txt> <background.csv|tsv|txt|NONE> <output_dir> [compound_target_table|NONE] [theme_keywords.csv|NONE] [p_cutoff] [padj_cutoff] [go_top_n] [kegg_top_n] [go_pool_per_ontology] [online_timeout_seconds] [input_set_label] [kegg_snapshot_dir|NONE] [kegg_fallback_mode]
```

Invoke the module once with `input_set_label="Drug-disease intersection targets"` and once with `input_set_label="Cytoscape core targets"`. The module writes ID-mapping tables, complete and filtered GO/KEGG tables, redundancy decisions, display tables, figures when eligible terms exist, optional compound-target links, a parameter table, session information, and QC.

If a downstream network or Sankey diagram uses a KEGG display table, explicitly record whether it came from the intersection-target analysis or core-target analysis. Do not merge both KEGG tables unless the user requests a separate comparison analysis.
