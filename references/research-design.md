# Research design and stage control

## Study definition

Define the research question before acquisition. Record whether the study concerns a full formula, selected herbs, a component group, a disease, a disease subtype, or a symptom concept. Do not broaden or narrow the disease without saying so.

Separate three kinds of decisions:

- **User decisions:** explicitly requested sources, thresholds, algorithms, figures, and exclusions.
- **Proposed defaults:** conventional or practical settings offered for confirmation.
- **Technical inferences:** format, encoding, column mapping, or harmless implementation choices made to continue work.

Store the decision type and rationale in the run record when it materially affects target selection or interpretation.

## Snapshot model

Use three data states:

1. `raw`: unchanged source exports or API responses.
2. `curated`: normalized and screened records with inclusion status and reasons.
3. `analysis`: frozen target sets and derived tables used for statistics and figures.

Never regenerate a figure from a live database call. A new database query creates a new raw snapshot and, after review, a new analysis snapshot.

## Stage status

Each stage should report one of:

- `ready`: required inputs are valid.
- `needs_review`: plausible results exist but a scientific choice is unresolved.
- `blocked`: the stage cannot produce a valid result.
- `complete`: outputs and quality-control checks are recorded.

## Expected project layout

```text
project/
|-- project-config.yaml
|-- 00_intake/
|-- 01_raw/
|-- 02_curated/
|-- 03_analysis/
|-- 04_figures/
|-- 05_report/
`-- logs/
```

Use descriptive stable filenames rather than embedding every parameter into the name. Store exact parameters in the configuration and manifest.

## Safe continuation

Before continuing from an existing project, identify the latest valid snapshot, verify expected columns and non-empty key sets, and state which stages will be reused. Do not overwrite user-supplied files. Derived outputs may be versioned by run identifier or date.

