# Two-pass STRING and Cytoscape core-target protocol

## First STRING network

Use the frozen drug-disease intersection gene list as the submitted protein list.

- Database: STRING.
- Species: `Homo sapiens`, taxon `9606`.
- Minimum required interaction score: `0.4`, corresponding to API `required_score = 400`.
- Network: interactions among submitted targets; do not add new neighbor proteins unless the user explicitly requests expansion.
- Hide disconnected nodes in the displayed network.
- Preserve the submitted target list, STRING mapping result, unmatched identifiers, access date, STRING version, and query parameters.

Export and retain at least:

1. high-resolution STRING network PNG;
2. short tabular interaction TSV suitable for Excel and Cytoscape.

Task-specific STRING download URLs contain temporary identifiers and must not be stored as reusable skill URLs. Store the downloaded files and the reproducible query parameters instead.

## Cytoscape import and analysis

Import the first-pass short TSV into Cytoscape as an undirected network. Map the protein A and protein B name or identifier columns to source and target. Preserve the STRING combined score as an edge attribute.

Run Cytoscape NetworkAnalyzer as an undirected network and export the node table with at least:

- node name or standardized gene symbol;
- Betweenness Centrality;
- Closeness Centrality;
- Degree.

Record the Cytoscape version, NetworkAnalyzer mode, node-name field, and whether the imported network contains multiple connected components. Do not calculate thresholds on rows missing one of the three metrics.

## Core-target rule

Calculate the median separately for Betweenness, Closeness, and Degree over valid analyzed nodes. A node is a core target only when all three conditions are true:

```text
Betweenness > median(Betweenness)
AND Closeness > median(Closeness)
AND Degree > median(Degree)
```

The comparison is strict `>`. Report the three exact median values, valid node count, excluded missing-metric rows, and final core-target count. Save the full node decision table as well as the unique core-target list.

Use [scripts/04_filter_core_targets_from_cytoscape.R](../scripts/04_filter_core_targets_from_cytoscape.R) when a Cytoscape node-table CSV or TSV is available.

## Second STRING network

Submit the resulting core-target gene list to STRING and repeat the first-pass settings:

- `Homo sapiens (9606)`;
- minimum interaction score `0.4` / `required_score = 400`;
- no added neighbor proteins;
- hide disconnected nodes;
- export high-resolution PNG and short interaction TSV.

Keep the first and second STRING outputs distinct. The second network is the core-target PPI network and must not overwrite the complete intersection-target PPI network.

Suggested filenames:

- `01_STRING_intersection_hires.png`
- `01_STRING_intersection_short.tsv`
- `02_Cytoscape_node_metrics.csv`
- `03_core_target_thresholds.csv`
- `04_core_targets.csv`
- `05_STRING_core_hires.png`
- `05_STRING_core_short.tsv`

