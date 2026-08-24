# Stage 02 and 03 schemas

## Stage 02 normalized source input

Each database export should be converted to CSV or TSV and contain `Gene_symbol`. Preserve all source-native columns. Add the following canonical columns when available:

| Column | Meaning |
|---|---|
| `Gene_symbol` | UniProt- or HGNC-standardized human gene symbol |
| `Score` | Source-native association score |
| `Evidence_type` | Curated, direct, inferred, genetic, literature, therapeutic, or other evidence |
| `Is_curated` | Explicit true or false when the source provides this distinction |
| `Association_type` | Direct, indirect, inferred, action target, and similar source relationship |
| `Disease_id` | Confirmed database disease identifier |
| `Disease_name` | Database preferred disease name |
| `Source_record_id` | Stable source record identifier |
| `PMID` | Supporting publication identifier when available |
| `Access_date` | Retrieval date |

Source adapters may map native score columns such as `Relevance_score`, `Inference_score`, `Overall_association_score`, or `GDA_score` into `Score`, while retaining the original field.

## Stage 02 outputs

- `01_疾病靶点原始汇总.csv`: all normalized records from all sources with source-native fields retained.
- `02_疾病靶点筛选决策.csv`: one row per source record with source rule, cutoff, pass status, and reason.
- `03_疾病靶点来源追踪.csv`: one row per gene and database after screening, including score and evidence summaries.
- `04_疾病靶点候选全集.csv`: unique broad candidate targets passing source-specific screening.
- `05_疾病靶点推荐集.csv`: curated or direct targets plus targets supported by at least two databases.
- `06_疾病靶点筛选统计.csv`: raw and retained counts, scored records, median or other cutoff, and missing-field warnings by source.
- `07_疾病靶点数据库统计.png` and `.pdf`: raw and retained unique-target counts by source.

## Stage 03 inputs

- a frozen drug compound-target table or unique drug-target table;
- either the broad candidate or recommended disease-target set;
- the detailed disease source-tracking table;
- study drug or formula label and confirmed disease label.

## Stage 03 outputs

- `01_药物疾病集合成员表.csv`: union of drug and disease targets with membership flags.
- `02_药物疾病交集靶点.csv`: unique intersecting targets.
- `03_交集靶点来源追踪表.csv`: intersecting targets with herbs, compounds, drug sources, disease databases, evidence tier, scores, and record counts.
- `04_药物疾病靶点韦恩图.png` and `.pdf`.
- `05_交集统计.csv`: drug, disease, union, and intersection counts plus the disease-set choice.

Use CSV content with comma delimiters when the extension is `.csv`. Do not write tab-delimited content under a `.csv` filename.

