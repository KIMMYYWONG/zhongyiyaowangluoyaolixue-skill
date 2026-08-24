# User intake

Use this guide at the start of a new study. Do not block useful planning merely because optional inputs are missing.

## Minimum inputs to begin

Ask for or infer from provided materials:

1. Formula name or study label.
2. Complete herb list, preferably with dosage, processed form, and source prescription.
3. Disease name and the intended disease scope or subtype.
4. Intended output: exploratory analysis, thesis or manuscript figures, methods section, target selection for docking, or another downstream use.

The framework can be created before these are finalized. Data acquisition should wait when herb identity or disease identity is ambiguous.

## Strongly useful materials

- Existing compound-target tables or database exports.
- Existing disease-target exports from GeneCards, OMIM, DisGeNET, Open Targets, or other specified sources.
- Screenshots or examples of preferred figure styles.
- Journal image specifications, language, dimensions, file formats, and resolution.
- Previous project folders whose table structures or naming conventions should be preserved.
- Rules previously agreed with the user, such as thresholds, evidence sources, manual additions, exclusions, and core-target algorithms.

## Intake record

Create a short intake record containing:

| Field | Required before acquisition | Notes |
|---|---:|---|
| Study label | Yes | Stable human-readable name |
| Formula and herbs | Yes | Include processing form and dosage if relevant |
| Disease query | Yes | Preserve the user's original wording |
| Confirmed disease ID | Yes | Confirm after showing plausible candidates |
| Species | Yes | Default proposal may be human, but record confirmation |
| Data sources | Yes | Distinguish requested sources from available sources |
| Screening rules | Yes | Mark user-specified versus proposed defaults |
| Intended use | Recommended | Controls output depth and figure standards |
| Figure language/style | Before final figures | Can remain undecided during analysis |

## Default proposals

If the user has not specified a value, present defaults as proposals rather than silently applying them. Common exploratory proposals include human targets, `OB >= 30`, `DL >= 0.18`, BH-adjusted enrichment, and high-confidence STRING interactions. These values must remain editable and must be written into the project configuration if used.

