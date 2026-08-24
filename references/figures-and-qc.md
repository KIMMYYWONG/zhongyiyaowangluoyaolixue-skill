# Figures and quality control

## Figure-data separation

Every figure must name its frozen input table and figure parameters. Style-only changes must reuse the same analysis snapshot. If filtering, Top N, or node selection changes, record it as a figure-selection change; if the underlying target set changes, create a new analysis snapshot.

## Core figure modules

Support independent generation of:

- drug-disease overlap visualization;
- herb-compound-target network;
- PPI network;
- core-target ranking plot;
- GO BP, CC, and MF plots, separate or combined;
- pathway enrichment plot;
- optional source-overlap, mapping-loss, or evidence-summary plots when they improve auditability.

Do not force every project to produce every figure. Choose figures that answer the research question and remain readable at the requested size.

## Publication outputs

For final figures, record language, dimensions, font, color palette, Top N, ordering, label wrapping, resolution, and output format. Prefer both a preview raster and an editable or vector format when supported. Always save the exact plotted data.

## Quality-control checks

Before delivery, verify:

- labels and titles match the confirmed formula and disease;
- legends, colors, node shapes, and abbreviations are defined;
- counts in figures, tables, and reports agree;
- no labels are clipped or unreadable at final size;
- Top N and display subsets are disclosed;
- p-values and adjusted p-values are not mislabeled;
- output filenames listed in the report actually exist;
- actual runtime parameters, rather than default text, appear in methods and logs.

Visually inspect generated figures. Numerical success alone is not sufficient quality assurance.

