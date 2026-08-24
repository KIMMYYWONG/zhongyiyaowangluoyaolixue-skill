---
name: network-pharmacology-research
description: Conduct staged, auditable network pharmacology research for traditional Chinese medicine formulas or herbs and diseases, including input normalization, drug and disease target selection, intersection analysis, PPI, enrichment, figures, and reproducible reporting. Use when the user wants to design, run, inspect, revise, or visualize a network pharmacology study; do not use for molecular docking alone or for literature summaries without data analysis.
---

# Network Pharmacology Research

Build a traceable research record, not only a set of figures. Preserve source data, screening decisions, identifier mappings, parameters, exclusions, and figure inputs so that every conclusion can be audited and rerun.

## Choose the working mode

- **Research design:** define the disease entity, formula or herbs, species, sources, thresholds, outputs, and intended publication use.
- **Staged analysis:** run one stage at a time and pause at the required confirmation gates.
- **Imported-data analysis:** accept user-supplied compound, target, disease, PPI, or enrichment tables and continue from the latest valid stage.
- **Figure-only revision:** redraw from a frozen analysis snapshot without querying databases or changing the selected target set.
- **Legacy quick run:** use the older one-command workflow only when the user explicitly prioritizes a rapid exploratory result over staged review.

Read [references/user-intake.md](references/user-intake.md) when starting a new study or identifying missing inputs. Read [references/research-design.md](references/research-design.md) before fixing the analysis plan.

## Required stage gates

1. **Confirm the research object.** Normalize the formula, herb names, disease name, disease ontology identifier, species, and study scope.
2. **Review drug-side selection.** Report unmatched herbs, retained and excluded compounds, target provenance, mapping losses, and the proposed final drug-target set.
3. **Review disease-side selection.** Show candidate disease entities before choosing one. Report each database separately, its evidence fields, normalization losses, deduplication, and the proposed final disease-target set.
4. **Freeze the target snapshot.** After approval or a clearly stated user instruction to proceed, calculate intersections, PPI, topology, and enrichment from immutable versioned tables.
5. **Freeze the analysis snapshot.** Draw figures only from saved tables and saved parameters. A style-only rerun must not re-fetch or silently alter data.
6. **Deliver with provenance.** Include tables, figures, configuration, database access dates or release information when available, run log, quality-control summary, and methods-ready text.

Do not silently pass a gate when disease resolution is ambiguous, a requested herb is unmatched, identifiers fail extensively, the intersection is empty, or a downstream method would be invalid. Present the evidence and the smallest decision the user needs to make.

## Analysis rules

- Treat familiar thresholds such as `OB >= 30` and `DL >= 0.18` as configurable study choices, not universal truths.
- Keep every raw value needed to reproduce inclusion and exclusion decisions. Never replace the raw table with only the filtered table.
- Before target normalization, resolve compounds from different databases to a stable canonical chemical key in R, deduplicate a compound dictionary, assign the reviewed preferred name from the Dictionary of Natural Products, and preserve every source identifier. Stop on name-only identities or identifier/name conflicts rather than merging them silently.
- Standardize targets to human gene symbols plus stable identifiers when available. Record aliases, unmapped rows, duplicates, species, and mapping source.
- Never accept the first disease search result automatically when multiple plausible entities exist.
- Keep disease databases separate until their evidence fields have been normalized. Record whether the final rule is union, intersection, multi-source support, or a score-based rule.
- When a requested record cannot be obtained through a structured export or ordinary lookup, use browser or computer webpage control to inspect the requested database interface. Preserve the same source, query, species, thresholds, and inclusion rules. Do not substitute a database, expand synonyms, relax a cutoff, or alter the frozen protocol without the user's explicit approval.
- Distinguish PPI display nodes, topology-ranked candidates, evidence-supported priority targets, and final experimental candidates.
- Record enrichment background, identifier conversion success, multiple-testing method, selection rule, ranking rule, and term-redundancy handling.
- Run GO and KEGG twice with `clusterProfiler`: one independent analysis for the frozen drug-disease intersection targets and one for the Cytoscape core targets. Both runs use Symbol-to-Entrez conversion, `enrichGO`/`enrichKEGG`, human organism settings, BH adjustment, and complete result capture before display filtering. Keep their directories and provenance separate; never pool their targets or results. Do not substitute Metascape unless the user explicitly changes this method.
- Avoid causal or clinical claims based only on network centrality or enrichment. Use language such as candidate target, associated pathway, or hypothesis for validation.

Read the stage-specific guidance only when that stage is active:

- Drug-side acquisition and screening: [references/drug-targets.md](references/drug-targets.md)
- TCMSP, SYMMAP, HERB 2.0, compound naming, and UniProt procedures: [references/drug-source-protocols.md](references/drug-source-protocols.md)
- Stage 01 table schemas and required quality-control fields: [references/stage-01-schema.md](references/stage-01-schema.md)
- Disease entity and target screening: [references/disease-targets.md](references/disease-targets.md)
- Disease database roles and source-specific screening: [references/disease-source-protocols.md](references/disease-source-protocols.md)
- Two-disease or phenotype comorbidity intersections: [references/comorbidity-intersection-protocol.md](references/comorbidity-intersection-protocol.md)
- Stage 02 and 03 disease and intersection schemas: [references/stage-02-03-schema.md](references/stage-02-03-schema.md)
- Intersection, PPI, core targets, and enrichment: [references/analysis.md](references/analysis.md)
- Two-pass STRING and Cytoscape core-target workflow: [references/ppi-core-protocol.md](references/ppi-core-protocol.md)
- GO/KEGG ID conversion, background selection, statistics, redundancy reduction, Top N figures, and pharmacology-theme summaries: [references/enrichment-protocol.md](references/enrichment-protocol.md)
- Browser/computer-assisted database screening and decision-freeze rules: [references/web-screening-and-change-control.md](references/web-screening-and-change-control.md)
- Formula-herb-ingredient-gene-disease Cytoscape tables, Excel ingestion with `openxlsx`, target-pathway Node/Edge tables, shared-ingredient numbering, and ingredient-gene-pathway Sankey diagrams: [references/network-sankey-protocol.md](references/network-sankey-protocol.md)
- Framework freeze, three-pass module validation, case-study replay, and GitHub release: [references/validation-and-release.md](references/validation-and-release.md)
- Figure generation and quality control: [references/figures-and-qc.md](references/figures-and-qc.md)
- Migration from the original pipeline: [references/legacy-migration.md](references/legacy-migration.md)

## Project configuration

For a new project, copy [assets/project-config.template.yaml](assets/project-config.template.yaml) into the project folder as `project-config.yaml`. Record only parameters actually used. If a value changes, create a new analysis snapshot rather than overwriting the previous provenance.

For disease-target aggregation, copy [assets/disease-source-policy.template.csv](assets/disease-source-policy.template.csv) into the project folder and adjust only when the study protocol requires a different source-specific rule.

For enrichment interpretation, copy [assets/pharmacology-theme-keywords.template.csv](assets/pharmacology-theme-keywords.template.csv) and revise its disease-specific keyword groups before running the module. Keywords are used to retrieve enriched terms for review, never to force a pathway into the conclusion.

## Minimum completion standard

A complete study must make it possible to answer:

- Which exact formula, herbs, disease entity, species, databases, dates, and thresholds were used?
- Which records were retained, excluded, unmatched, duplicated, or unmapped, and why?
- Which frozen tables generated each analysis and figure?
- Which choices were user-specified, which were proposed defaults, and which were inferred?
- Can another run reproduce the counts, target lists, statistical settings, and visual selections?

Do not start the real formula case study until the framework files, schemas, configurable decisions, and synthetic checks are complete. After framework freeze, replay the study from Stage 01 with the chosen formula and disease, validating each module three ways before proceeding to the next stage.
