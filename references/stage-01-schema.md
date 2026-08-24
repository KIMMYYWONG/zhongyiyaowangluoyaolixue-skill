# Stage 01 schemas: herbs, compounds, and targets

Use flat tables with one record per row. Do not rely on colors or worksheet names as the only source marker.

## `herbs.csv`

| Column | Meaning |
|---|---|
| `herb_input_cn` | User-provided Chinese herb name |
| `herb_normalized_cn` | Confirmed normalized Chinese name |
| `herb_pinyin` | Stable pinyin or database display name |
| `processed_form` | Processing form when relevant |
| `tcmsp_herb_id` | TCMSP herb identifier when available |
| `symmap_smhb_id` | SYMMAP SMHB identifier |
| `herb2_herb_id` | HERB 2.0 identifier |
| `match_status` | exact, ambiguous, unmatched, or pending |
| `review_note` | Resolution evidence or user decision |

## `compound_raw.csv`

One row per source herb-ingredient record.

| Column | Meaning |
|---|---|
| `source` | TCMSP, SYMMAP, or HERB2 |
| `source_herb_id` | Herb identifier in that source |
| `source_ingredient_id` | MolID, SYMMAP ingredient ID, or HERB ID |
| `herb_normalized_cn` | Confirmed herb name |
| `herb_pinyin` | Pinyin or database display name |
| `compound_name_raw` | Compound name exactly as displayed |
| `ob` | Numeric OB value when available |
| `dl` | Numeric DL value when available |
| `pubchem_cid_raw` | Source PubChem CID |
| `cas_raw` | Source CAS value |
| `inchi_key_raw` | Source structure key when available |
| `evidence_type` | Database, reference, experiment, prediction, or other source class |
| `source_url` | Record or detail-page URL |
| `access_date` | Retrieval date |

## `compound_decisions.csv`

One row per reviewed source compound identity. Include raw identifiers plus:

| Column | Meaning |
|---|---|
| `source_compound_key` | Source-side key generated from InChIKey, PubChem CID, CAS, or temporary normalized name |
| `compound_key` | Reviewed canonical stable key shared by records representing the same chemical |
| `compound_name_preferred` | Preferred normalized name |
| `name_authority` | Dictionary of Natural Products evidence used for the preferred name |
| `pubchem_cid` | Accepted PubChem CID |
| `cas` | Accepted CAS identifier or delimited set when identity is verified |
| `inchi_key` | Accepted InChIKey |
| `identity_match_method` | Structure, CID, CAS, or reviewed name matching |
| `identity_match_status` | exact, ambiguous, unmatched, or pending |
| `screen_rule` | Source-specific rule applied |
| `include_compound` | true or false |
| `include_reason` | Threshold, supplementary evidence, manual addition, or exclusion reason |

`source_compound_key` must be unique in the decision table, but multiple source keys may resolve to the same canonical `compound_key`. Keys beginning with `NAME:` are temporary unresolved keys and cannot enter the frozen output. TCMSP MolID, SYMMAP ingredient ID, and HERB ID remain source identifiers rather than canonical chemical keys.

## `01a_compound_identity_dictionary.csv`

One row per canonical stable compound key after cross-database resolution. It records the DNP preferred name, accepted PubChem CID/CAS/InChIKey, identity method and status, conflict flag, contributing database count and list, source ingredient IDs, TCMSP MolIDs, herbs, and raw names. The accompanying `01b_compound_cross_source_dedup_audit.csv` keeps the multi-source or conflicting dictionary entries for review; `01_compound_screening_decisions.csv` keeps the complete one-row-per-source audit with source key, canonical key, eligibility, identity readiness, DNP naming readiness, and inclusion decision.

## `target_raw.csv`

One row per source ingredient-target relationship.

| Column | Meaning |
|---|---|
| `source` | Source database |
| `source_herb_id` | Source herb ID |
| `source_ingredient_id` | Source ingredient ID |
| `compound_key` | Normalized compound key if resolved |
| `target_name_raw` | Target text exactly as displayed |
| `target_id_raw` | Source target ID when available |
| `evidence_type` | Relationship evidence category |
| `evidence_reference` | Reference or experiment identifier |
| `source_url` | Target or ingredient detail URL |
| `access_date` | Retrieval date |

For HERB rows, also retain `pmid`, `doi`, or another traceable citation when available. A final HERB supplement row must be traceable to literature evidence.

## `target_uniprot_mapping.csv`

| Column | Meaning |
|---|---|
| `target_name_raw` | Original target name |
| `organism_requested` | Usually Homo sapiens |
| `uniprot_accession` | Accepted UniProt accession |
| `uniprot_entry_name` | UniProt entry name |
| `protein_name_recommended` | Recommended protein name |
| `gene_symbol_primary` | Accepted primary gene symbol |
| `gene_names_alternative` | Alternative symbols or names |
| `reviewed_status` | reviewed or unreviewed |
| `organism_result` | Organism of matched entry |
| `mapping_method` | ID mapping, exact name, alias, or reviewed manual mapping |
| `mapping_status` | exact, alias, ambiguous, unmatched, or non-human |
| `mapping_note` | Evidence and resolution notes |

## `compound_target_final.csv`

One row per retained normalized herb-compound-target relationship. Keep the user's legacy-compatible fields first:

1. `herb_pinyin`
2. `pubchem_cid`
3. `tcmsp_molid`
4. `compound_name_preferred`
5. `gene_symbol_primary`
6. `cas`
7. `source`

Then append `source_herb_id`, `source_ingredient_id`, `uniprot_accession`, `evidence_type`, `source_url`, `access_date`, `compound_decision`, and `target_mapping_status`.

For HERB relationships, also append `herb_selection_tier` with one of:

- `tcmsp_target_overlap`;
- `representative_ingredient_target`;
- `not_selected_volume_control`.

Store `representative_ingredient_reason` and `evidence_reference` so that the supplementary choice can be audited.

## `04a_compound_target_deduplicated.csv`

One row per unique `herb_normalized_cn + compound_key + gene_symbol_primary` biological relationship. Collapse contributing databases, source herb IDs, source ingredient IDs, references, and URLs into delimited provenance fields. Keep `compound_target_final.csv` as the source-preserving compatibility export; use the deduplicated table for downstream unique-relation analysis.

## Required summaries

Produce per-herb and per-source counts for raw compounds, eligible compounds, normalized compounds, raw targets, accepted human targets, ambiguous mappings, unmapped targets, and final unique targets. Preserve shared compounds and shared targets rather than assigning them to only one herb.

Use [scripts/01_drug_target_normalize_filter.R](../scripts/01_drug_target_normalize_filter.R) after the source-specific raw tables, compound-name decisions, and UniProt mapping review are available:

```text
Rscript 01_drug_target_normalize_filter.R <compound_raw.csv|tsv> <target_raw.csv|tsv> <compound_decisions.csv|tsv|NONE> <uniprot_mapping.csv|tsv> <output_dir> [ob_min=30] [dl_min=0.18]
```

The script re-evaluates source-specific rules instead of trusting prefiltered input, resolves stable source keys to canonical compound keys, blocks identity conflicts and name-only keys, requires a Dictionary of Natural Products preferred name, accepts only reviewed human UniProt mappings, handles HERB literature and volume-control tiers, and writes the compound dictionary, cross-source audit, source-preserving final table, and deduplicated biological relationship table.
