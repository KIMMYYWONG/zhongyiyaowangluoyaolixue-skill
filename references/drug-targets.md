# Drug-side compounds and targets

## Herb normalization

Preserve the original herb text and map it to a normalized name. When available, record Latin or pharmacopoeial name, processed form, medicinal part, dosage, and source prescription. Report every unmatched or multiply matched herb before acquisition continues.

Do not silently treat processed and unprocessed forms, similarly named species, or substitute herbs as identical.

## Compound acquisition

For every compound record, preserve:

- source database and record identifier;
- herb association;
- compound name and stable identifier;
- raw screening fields such as OB and DL;
- acquisition date or database release when available;
- inclusion status and reason;
- whether the compound was retained by threshold, literature evidence, pharmacopoeial relevance, or an explicit user rule.

If the user asks to add known active compounds that fail conventional thresholds, keep them in a separately identified evidence class rather than altering their raw OB or DL values.

For the user's established workflow, treat TCMSP as the primary screening source and SYMMAP plus HERB 2.0 as supplementary sources. Read [drug-source-protocols.md](drug-source-protocols.md) before acquiring or merging these records. Do not apply the TCMSP `OB >= 30` and `DL >= 0.18` pair silently to another database when the user has specified a different source-specific rule.

## Target normalization

Before joining targets, resolve cross-database compound records to canonical stable chemical keys in R and assign the reviewed Dictionary of Natural Products preferred name. Then preserve each raw target name and map it through UniProt to a human primary gene symbol plus stable accession. Record species, mapping source, alias resolution, mapping status, and reasons for exclusion.

Produce at least:

- raw compound-target table;
- curated compound-target decision table;
- mapping exception table;
- deduplicated drug-target table;
- herb-level counts before and after screening.

Use the schemas in [stage-01-schema.md](stage-01-schema.md). Retain the user's legacy six-column export order for compatibility, but add a source column and keep OB, DL, identifiers, mapping evidence, inclusion status, and reasons in the audit tables.

## Drug-side review gate

Before freezing the drug-target set, show:

- herbs requested, matched, unmatched, and ambiguous;
- compound counts before and after each rule;
- target counts before normalization, after normalization, and after deduplication;
- compounds or herbs contributing no retained targets;
- any manual additions or exclusions and their evidence class.

Do not continue silently if a requested herb is absent, the mapping loss is substantial, or the final target set is empty.
