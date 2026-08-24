# Intersection, PPI, core targets, and enrichment

## Intersection

Calculate the intersection only from frozen normalized drug and disease target tables. Save the membership table, not only a Venn image. Include source support, drug-side connections, disease evidence, and mapping identifiers when available.

Use the schemas in [stage-02-03-schema.md](stage-02-03-schema.md). The intersection output must identify the selected disease-target set, drug databases or herbs, disease databases, evidence tier, and source-screening decision behind every overlapping gene.

If the intersection is empty or unexpectedly small, stop and audit disease identity, species, identifier mapping, thresholds, and source coverage before downstream analysis.

## PPI

Record species, STRING or other network source, database version or access date, interaction score threshold, network type, and treatment of isolated nodes and duplicated edges. Save submitted identifiers, returned mappings, excluded identifiers, raw interactions, and cleaned edges.

Do not describe a display-only subset as the complete PPI network.

For the user's established workflow, use the two-pass STRING and Cytoscape procedure in [ppi-core-protocol.md](ppi-core-protocol.md).

## Core-target selection

Keep all computed topology metrics. Distinguish:

- all PPI nodes;
- the subset displayed in a figure;
- topology-ranked candidate targets;
- targets supported by multiple evidence dimensions;
- final targets selected for docking or experiments.

Record the ranking or filtering algorithm. Degree alone may be used for an exploratory ranking, but should not automatically become a biological conclusion. When combining metrics, document scaling, weights, tie handling, and cutoff.

The user's default Cytoscape rule is strict simultaneous filtering: Betweenness greater than its median, Closeness greater than its median, and Degree greater than its median. Do not replace strict `>` with `>=` or with an OR rule.

## GO and pathway enrichment

Use `clusterProfiler` for both GO and KEGG enrichment. Run the workflow independently for (1) the frozen drug-disease intersection targets and (2) the Cytoscape core targets. Convert each set to Entrez IDs separately with `org.Hs.eg.db`, retain mapping losses, and save complete results before filtering or plotting. Never combine the two target lists or result tables.

Record:

- analyzed target set and gene count;
- submitted gene identifiers, rejected identifiers, and accepted Entrez count;
- organism and annotation database;
- clusterProfiler background universe and statistical settings;
- ontology or pathway database and release when available;
- statistical test, p-value cutoff, adjustment method, adjusted cutoff, and q-value rule;
- ranking rule used to select displayed terms;
- term-redundancy reduction method, if any.

Save the complete enrichment result before selecting terms for figures. A figure's Top N is a display choice and must not redefine statistical significance.

If enrichment fails or returns no eligible terms, preserve the diagnostic evidence and do not fabricate an empty-looking figure.

Use the full procedure in [enrichment-protocol.md](enrichment-protocol.md). For the user's default publication workflow, analyze the Cytoscape-selected core targets, display 15 GO terms in total using balanced BP/MF/CC selection when possible, and display 30 nonredundant KEGG pathways. Keep raw `P < 0.05` eligibility because it was explicitly requested, but always report BH-adjusted values and identify which displayed terms also satisfy adjusted `P < 0.05`.
