# Framework validation, case replay, and release

## Work order

Use this fixed order:

1. finish and freeze the reusable framework;
2. run synthetic positive and negative fixtures without using the real study as a design shortcut;
3. start the user-confirmed case study from Stage 01 only after framework freeze;
4. complete one module, perform all three checks, and freeze its outputs before the next module;
5. revise the reusable framework whenever the case exposes a general defect, then rerun affected downstream modules;
6. publish to the user-confirmed GitHub repository only after the complete case replay passes.

The case study is a validation case, not hard-coded skill content. Formula-specific herb names, abbreviations, disease terms, thresholds, and conclusions belong in the project configuration and case output, not in reusable scripts.

## Three checks required for every module

### Check 1: schema and scientific-rule review

Confirm required columns, identifiers, species, thresholds, inclusion logic, exclusion reasons, source provenance, and user-selected interpretation. Review row-level decision tables rather than only counts or figures.

### Check 2: automatic integrity and negative tests

Run the module's automatic QC and at least one deliberately invalid fixture. Confirm that missing identifiers, invalid enrichment backgrounds, mixed shared-ingredient modes, blank endpoints, duplicate edges, or empty intersections stop with an informative error rather than generating plausible-looking output.

### Check 3: cross-module and reproducibility review

Rerun from the same frozen inputs and compare counts, stable identifiers, hashes where available, and figure-input tables. Verify that downstream genes are subsets of the correct upstream set, displayed Top N rows exactly match saved display tables, and style-only reruns do not change data.

Record every pass, failure, correction, rerun, reviewer, and timestamp in [module-validation-log.template.csv](../assets/module-validation-log.template.csv). A module advances only when all three checks pass.

## Framework-freeze checklist

- Every stage has declared inputs, outputs, configurable decisions, stop conditions, and QC fields.
- No script contains a case-specific absolute working directory or fixed formula, disease, herb list, or conclusion.
- Complete/raw tables are preserved separately from filtered and display tables.
- The `complete` versus `filtered` shared-ingredient decision is explicit and immutable within a snapshot.
- Online failures are recorded without silently replacing data sources.
- Example fixtures contain synthetic or generic data only.
- Required packages and optional plotting packages are documented; absent optional packages do not corrupt tables.

## GitHub release gate

Before publishing, confirm the target repository, branch, visibility, license, and whether example outputs may be public. Exclude credentials, cookies, temporary download URLs, local absolute paths, private literature files, and raw database exports whose terms do not permit redistribution.

Commit the reusable skill, templates, scripts, generic fixtures, methods documentation, and validation log. Keep real case data in a separate example directory only when the user approves its publication. Tag the release only after a clean checkout can validate the skill and reproduce the permitted example workflow.
