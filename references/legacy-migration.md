# Migration from the original pipeline

Treat the original `network-pharmacology-pipeline` as a useful exploratory implementation. Preserve it unchanged while building the staged workflow.

## Reusable components

- TCMNP or TCMSP cache acquisition.
- Initial herb-compound-target extraction.
- Open Targets API access.
- STRING API access.
- Basic intersection, PPI metrics, GO, KEGG, and plotting code.
- One-command orchestration as an optional legacy quick mode.

## Changes required before staged use

- Replace environment-specific executable paths with detected or configured paths.
- Preserve Chinese study names with a stable project identifier instead of reducing unsupported names to `result`.
- Add explicit disease-candidate review before fetching final disease targets.
- Move OB, DL, disease score, STRING score, Top N, and figure settings into project configuration.
- Make reports read actual runtime settings rather than fixed prose.
- Align report file references with files actually generated.
- Add validation for unmatched herbs, empty target sets, empty intersections, empty PPI results, mapping losses, and empty enrichment results.
- Separate acquisition, curation, analysis, and figure commands so a later stage can be rerun independently.
- Handle compounds shared by multiple herbs without assigning them silently to only one herb for coloring.
- Preserve raw API responses or source exports and record access dates.

## Migration sequence

1. Define standard schemas and project configuration.
2. Wrap existing acquisition scripts so they write raw and curated outputs separately.
3. Add review summaries and confirmation gates.
4. Refactor intersection, PPI, and enrichment to consume frozen standardized tables.
5. Split each figure into an independent renderer.
6. Generate manifest, QC summary, methods text, and output index from actual files and parameters.
7. Retain the legacy one-command runner only as a clearly labeled exploratory mode.

For legacy compound-target workbooks with one herb per worksheet and the first three columns representing ingredient ID, ingredient name, and gene symbol, use [scripts/00_legacy_multisheet_compound_target_adapter.R](../scripts/00_legacy_multisheet_compound_target_adapter.R) with an explicit sheet-index-to-herb map. The explicit map prevents corrupted or locale-dependent worksheet labels from silently changing herb identity. Treat the resulting table as an imported legacy snapshot; it does not supply missing OB/DL, stable chemical identity, DNP naming, or UniProt evidence for a formal Stage 01 freeze.
