# Release checklist

- [x] `SKILL.md` frontmatter and body validated with equivalent checks.
- [x] `agents/openai.yaml` parsed and required interface fields verified.
- [x] All R scripts parse successfully under R 4.6.1.
- [x] Main KOA–Sarcopenia comorbidity replay completed.
- [x] Positive, negative, and deterministic repeat checks recorded.
- [x] Real database exports and generated case outputs excluded by default.
- [x] Target repository: `KIMMYYWONG/network-pharmacology-research-skill`.
- [x] Repository visibility: private.
- [x] Target branch: `main`; repository was empty before the initial release.
- [x] License: private/no-license use for this initial release.
- [x] Real case data and database exports are excluded; only source, templates, synthetic examples, and selected human-readable validation reports are published.
- [x] Connector publication scope reviewed; initial `main` commit prepared and pushed.

The official Python `quick_validate.py` could not run in the current bundled Python because PyYAML is absent. The same validations were run with R `yaml`: allowed frontmatter keys, required fields, name format/length, description length/characters, and unfinished TODO checks.
