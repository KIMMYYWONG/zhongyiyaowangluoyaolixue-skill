# Comorbidity disease-target intersection

Use this protocol when the study concerns two explicitly confirmed diseases or phenotypes and the intended comorbidity target set is their intersection.

## Required order

1. Resolve and confirm each disease independently, including its database identifiers and species.
2. Run source-specific screening independently for disease A and disease B. Never pool their raw database scores.
3. Freeze the selected target tier for each disease. Both diseases must use the same declared tier, such as `broad_candidate` or `evidence_supported`, unless the protocol explicitly justifies an asymmetric choice.
4. Intersect standardized human gene symbols and save the union membership table, common targets, disease-specific source tracking, counts, and QC.
5. Use the frozen comorbidity intersection as the disease-side input for the later drug-disease intersection.

Do not call a union of two disease lists a comorbidity intersection. Do not mix one disease's raw list with the other's filtered list. If the intersection is empty or unexpectedly small, stop and audit entity identity, species, source availability, screening tier, and identifier mapping.

Run [scripts/02b_disease_comorbidity_intersection.R](../scripts/02b_disease_comorbidity_intersection.R):

```text
Rscript 02b_disease_comorbidity_intersection.R <disease_a_targets> <disease_b_targets> <disease_a_tracking|NONE> <disease_b_tracking|NONE> <output_dir> <disease_a_label> <disease_b_label> [set_label]
```

The overlap figure is descriptive. The saved membership and source-tracking tables, not the figure, are the frozen analytical record.

