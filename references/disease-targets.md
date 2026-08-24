# Disease entity and disease targets

## Disease resolution

Preserve the user's disease wording. Search the selected ontology or database and present plausible candidates with identifiers, preferred names, synonyms, and short distinctions. Require confirmation when more than one candidate could materially change the study.

Never select the first search hit solely because it ranks first.

## Source handling

Keep each disease database in its native evidence model before combining:

- source database and disease identifier;
- source target identifier and normalized gene identifier;
- original score, evidence type, or provenance field;
- species and target biotype when available;
- mapping status;
- inclusion status and reason;
- acquisition date or database release when available.

Do not compare unlike scores as though they share a common scale. If normalized scores are created, document the transformation and retain the originals.

Read [disease-source-protocols.md](disease-source-protocols.md) before screening GeneCards, DrugBank, CTD, Open Targets, DISGENET, and OMIM exports. Use the source-specific policy rather than applying one pooled median to scores from different databases.

## Combination strategies

Record the exact final rule. Supported research choices may include:

- union of all eligible source targets;
- intersection of specified sources;
- support from at least a chosen number of sources;
- source-specific thresholds followed by union;
- evidence-weighted ranking followed by a declared cutoff;
- a user-supplied curated target list.

The skill should not imply that one strategy is universally correct. The choice should match the protocol and be visible in the methods report.

The recommended default produces two frozen sets:

- a **broad candidate set** containing every target that passes its own source rule;
- an **evidence-supported recommended set** containing targets with curated or direct evidence, or support from at least two independent databases.

Keep both sets and state which one is used downstream. This allows conventional median-based screening without discarding stronger evidence that has no comparable numeric score.

## Disease-side review gate

Before freezing the disease-target set, show:

- confirmed disease entity and identifier;
- per-source raw, eligible, normalized, and unique target counts;
- overlaps among sources when multiple sources are used;
- unmapped and excluded records;
- the proposed combination rule and resulting count;
- warnings about unavailable requested databases or incomplete evidence.
