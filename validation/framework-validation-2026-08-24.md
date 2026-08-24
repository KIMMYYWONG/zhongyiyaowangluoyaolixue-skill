# Framework validation — 2026-08-24

## Scope

This validation covers only the reusable network-pharmacology skill framework. No real formula, herb list, disease dataset, or study conclusion was used.

## Result

Status: **PASS — framework freeze candidate**.

All eight runnable modules passed three checks: schema/scientific-rule review, automatic integrity with an expected-failure test, and frozen-input reproducibility or downstream-subset review.

## Check 1 — structure and scientific rules

- `SKILL.md`, UI metadata, 16 routed references, five templates, eight scripts, and generic fixtures are present.
- All local Markdown links resolve.
- All eight R scripts parse successfully.
- No case-specific formula, herb list, desktop path, or unfinished placeholder occurs in the reusable skill.
- The skill frontmatter passed an equivalent of `quick_validate.py`. The bundled Python validator itself could not start because its environment lacks PyYAML; this is a validator-environment limitation, not a skill-frontmatter failure.
- Stage rules match the framework configuration: source-specific drug screening, source-specific disease screening, frozen intersection, two-pass PPI, strict three-metric core rule, auditable GO/KEGG, explicit shared-ingredient mode, and frozen-table Sankey construction.

## Corrections made during validation

1. Disease-source input now stops when its source policy is missing or ambiguous, when a screening mode is unsupported, or when policy parameters are invalid.
2. An empty drug–disease intersection now stops before Venn, PPI, or enrichment.
3. A Cytoscape result with no targets satisfying the strict AND rule now saves its decision tables and stops before the second STRING run.
4. Sankey pathway records now preserve the frozen KEGG p-value and calculated `-log10(P)` in both triple and node tables.
5. Browser/computer-assisted screening is now the required fallback when a requested database can only be completed through its webpage; access-method fallback cannot change the source or scientific rules.
6. Both GO and KEGG use `clusterProfiler`, matching the confirmed senior workflow; the intersection-target and core-target sets are analyzed independently and Metascape is not an automatic substitute.
7. Network, target-pathway, and Sankey modules now accept explicitly selected Excel worksheets through `openxlsx::read.xlsx`; target-pathway output preserves the senior `from_node`/`to_node` and `node`/`type` sheets while adding canonical edge, node, QC, and provenance tables. Core and intersection KEGG results are forced to remain separate.
8. The Sankey module now attaches `ggsankeyfier` before using its registered `sankeynode` statistic; an actual 600-dpi PNG and vector PDF generation test passed.
9. Stage 01 now resolves source compounds to stable canonical chemical keys before joining targets, requires Dictionary of Natural Products preferred-name evidence, preserves all database IDs, blocks identity/key conflicts, standardizes targets through reviewed human UniProt mappings, and writes both source-preserving and cross-source-deduplicated relationship tables.
10. The two-herb case replay exposed and fixed partial-decision `NA` indexing, multi-ID strings being treated as unique keys, and insufficiently explicit QC for missing TCMSP/SYMMAP screening fields and HERB literature references. The standard positive fixture still passes after these corrections.

## Check 2 — positive and negative synthetic tests

Positive synthetic workflow:

- Stage 01: 3 retained compound–target rows and 3 unique drug targets.
- Stage 01 cross-source identity fixture: TCMSP and SYMMAP records with the same PubChem CID resolved to one compound-gene relationship while retaining both sources and both source ingredient IDs.
- Stage 02: 9 raw disease records, 4 broad targets, and 3 recommended targets.
- Stage 03: 1 frozen drug–disease intersection target.
- Stage 04: 4 valid topology rows and 2 strict-AND core targets.
- Stage 05: independent intersection-target and core-target runs wrote distinct `05A` and `05B` directories with correct input-set labels. Symbol-to-Entrez mapping and `clusterProfiler::enrichGO` completed; KEGG used `clusterProfiler::enrichKEGG` and recorded the current SSL retrieval failure without changing sources.
- Stage 06: both `complete` and `filtered` shared-ingredient modes passed all network QC checks.
- Stage 07: 5 unique ingredient–gene–pathway triples were produced.
- Stage 08: 2 KEGG pathways, 4 selected genes, 5 pathway-gene edges, and 5 ingredient-gene-pathway triples were produced. CSV and explicitly selected multi-sheet XLSX inputs gave the same relationship counts.

Eight deliberately invalid tests all stopped with informative errors:

1. nonnumeric Stage 01 threshold;
2. disease input without a matching source policy;
3. empty drug–disease intersection;
4. no node satisfying all three core-target thresholds;
5. core target absent from a supplied custom enrichment background;
6. missing herb code in network construction;
7. no matched ingredient–gene–pathway triple.
8. multi-sheet Excel input without an explicit worksheet selection.

The extended Stage 01 naming test also stopped when PubChem was supplied as the preferred-name authority instead of the Dictionary of Natural Products.

## Check 3 — reproducibility and module linkage

The primary frozen CSV outputs from Stages 01–08 were regenerated and compared by SHA-256. All eight module hashes matched. The extended Stage 01 deduplicated compound-target table also matched across repeat runs. Drug targets, disease targets, intersection targets, core targets, clusterProfiler-derived GO tables, network edges, Sankey triples, and target-pathway tables were stable across reruns. KEGG reproducibility remains conditional on successful access to the same online annotation service.

## Environment-limited behavior verified

- `ggsankeyfier 0.1.8` is installed with `gridBezier 1.1-1` and `vwline 0.2-4`; the Sankey module generated both a 600-dpi PNG and vector PDF.
- The active R runtime now provides the user-requested `openxlsx 4.2.5.2`; Excel write/read round-trip validation passed.
- Windows R emitted locale-startup warnings, but UTF-8 output and all tested CSV contents remained valid.

## Gate decision

The reusable framework may now be frozen for the next phase. The next phase is a separate case replay beginning at Stage 01; it must not retroactively hard-code case-specific names or conclusions into this framework.
