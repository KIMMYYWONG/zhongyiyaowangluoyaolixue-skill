# Cytoscape network and Sankey protocol

## Freeze and normalize the shared data backbone

The Cytoscape network and Sankey diagram must use the same frozen normalized compound-target relationships. Standardize herb names, standardized compound names, and human gene symbols; trim whitespace; remove empty values; and remove exact duplicate records.

Determine compound identity in this order when available: InChIKey, PubChem CID, CAS, then normalized standardized compound name. Do not merge records solely because two unstandardized display names look similar. Preserve original row identifiers and every identifier used in identity resolution.

Match the frozen drug-disease intersection genes against every herb's `Gene_symbol`. Save the row-level retained table with herb, original ingredient identifier, standardized ingredient name, gene, source database, and source row ID. This is the first manual-check table.

## Shared-ingredient assignment is a required decision

Choose exactly one mode for a study snapshot:

- `complete`: connect a retained shared ingredient to every herb containing that same compound in the full normalized herb-compound source table, even when a particular herb has no retained intersection-target row for it;
- `filtered`: connect the ingredient only to herbs that still have a retained ingredient-intersection-gene record.

Never combine edges generated under different modes. Record the selected mode in every output and audit table. In `complete` mode, require a frozen full herb-compound membership table or explicitly confirm that the supplied compound-target table is the complete membership source.

## Ingredient numbering

Resolve shared status before assigning node IDs. Sort deterministically by compound identity key.

- An ingredient belonging to only one herb receives that herb's lowercase code plus a sequence, such as `bs1`, `bs2`.
- An ingredient belonging to more than one herb receives one shared node ID, such as `S1`, `S2`.
- The same compound may never receive multiple ingredient nodes within the same frozen snapshot.

Save a node dictionary containing ingredient node ID, identity key, standardized name, shared flag, herb count, herb list, assignment mode, and retained genes.

## Cytoscape edge and node tables

Build a three-column edge table:

| source | target | interaction |
|---|---|---|
| formula | herb | `contains_drug` |
| herb | ingredient node | `contains_ingredient` |
| ingredient node | intersection gene | `targets` |
| disease | intersection gene | `associated_with` |

Build a node table with `term`, `type`, `label`, and `herbs`. Allowed types are `Formula`, `Drug`, `Ingredient`, `Gene`, and `Disease`. `term` is the unique import key; `label` is the displayed name. Every node occurs exactly once.

Use [scripts/06_build_network_tables.R](../scripts/06_build_network_tables.R). Import its edge table into Cytoscape as a Network with `source` and `target`, then import the node table as a Node Table using `term` matched to network `shared name`. Map color, shape, and size from `type`.

CSV, TSV, TXT, XLSX, and XLSM inputs are accepted. Excel files are read with `openxlsx::read.xlsx(check.names = FALSE)`. A workbook with more than one worksheet must identify the sheet explicitly as `file.xlsx::SheetName`; never silently use the first sheet. Record the runtime `openxlsx` version in provenance. Version `4.2.5.2` may be recorded as the legacy reference version, but do not claim it was used unless the current runtime actually reports that version.

## Target-pathway network compatibility

The senior workflow generates a pathway-gene edge table (`from_node`, `to_node`) and a node-type table (`node`, `type`) for Cytoscape. Preserve those compatibility sheets, but also write canonical tables using `source`, `target`, `interaction` and `term`, `type`, `label`. Use [scripts/08_build_target_pathway_network.R](../scripts/08_build_target_pathway_network.R) to:

1. read the frozen target set and KEGG result table;
2. filter `pvalue < 0.05`, rank ascending by p-value, and retain up to 30 pathways unless the frozen protocol says otherwise;
3. expand the pathway gene list and use `dplyr::inner_join()` to retain only genes in the explicitly selected target set;
4. optionally join the normalized ingredient-gene table to produce ingredient-gene-pathway triples;
5. deduplicate relations and export Cytoscape edge, node, QC, and provenance sheets.

Set `target_set` to exactly `intersection` or `core`. The target list and KEGG result must come from the same independent enrichment run. A core-target network must use the core-target KEGG output; an intersection-target network must use the intersection-target KEGG output. The legacy example loaded the intersection KEGG result and then intersected it with core targets; retain this only as historical evidence, not as the reusable default.

## Required double audit

First audit the row-level number-ingredient-gene table. Second compare raw, normalized, retained, and distinct combination counts. The automatic QC must verify:

- no empty edge endpoints;
- no duplicate `source-target-interaction` edges;
- every edge endpoint occurs in the node table;
- every node key occurs once;
- no node key has conflicting types;
- every shared ingredient has at least two assigned herbs;
- every unique ingredient has exactly one assigned herb;
- herb-ingredient edges conform to the selected assignment mode.

Stop before Cytoscape import when a critical check fails.

## Sankey framework

The default Sankey follows the user's established R workflow but reads frozen framework outputs:

```text
Ingredient -> core/intersection Gene -> nonredundant KEGG pathway
```

Expand each KEGG pathway's saved `Gene_symbols`, intersect them with the frozen ingredient-gene table, and create one unique ingredient-gene-pathway triple. Apply display limits only after saving the full triple table. Default display limits are the 30 ingredients with the most retained genes and 30 genes with the broadest retained connections; these are figure readability choices, not biological significance thresholds.

Use [scripts/07_build_sankey.R](../scripts/07_build_sankey.R). It accepts the same Excel sheet syntax as the network scripts. It always writes the full and displayed triple tables plus node metadata. If `ggsankeyfier` is installed, it also exports 600-dpi PNG and vector PDF. If the package is absent, retain the tables and record a plotting dependency warning rather than altering the data.

Pathway labels and `-log10(P)` values must come from one explicitly selected frozen KEGG display table. Record whether the table belongs to the intersection-target enrichment or core-target enrichment; do not mix the two analyses inside one Sankey diagram. Do not independently reselect pathways inside the Sankey script. For a single-molecule study, the ingredient layer may contain one node; for a formula, preserve the shared-ingredient node dictionary and herb list.
