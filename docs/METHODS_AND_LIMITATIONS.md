# Methods, reproducibility and limitations

## Discovery cohort
Site-balanced TCGA-LUAD 15 tumors and 15 normals; 7 overlapping TSS groups.
DESeq2 on raw STAR unstranded counts, model `~ tss + sample_type`.
Genes with BH-adjusted p <0.05 and absolute log2 fold change >1
are labeled exploratory differential-expression hits.

## Functional interpretation
GO Biological Process 2023 and Reactome 2022, GSEApy 1.3.2,
1,000 permutations, seed 42. ORA uses the tested-gene universe;
preranked GSEA uses the DESeq2 Wald statistic. Overlapping pathway
terms are not independent hypotheses.

## External validation
GSE32863 paired microarray study; 57 matched pairs after excluding
unpaired patients 3023 and 3035. Paired limma fixed-effect model.
One probe per symbol chosen by across-sample variance.
A gene 'replicates' if external BH FDR <0.05 and effect direction matches.
The direction-agreement denominator is the 2,600 original TCGA-significant
genes with external mappings, not all 4,469 discovery hits.

The 20 pathway tests transfer **TCGA-selected leading-edge gene sets**
into the external cohort; each score averages per-gene z-normalized
external expression, then applies a paired Wilcoxon test and BH
correction. These are **not** independent Reactome/GO GSEA re-estimates.
Leading-edge overlap and discovery-based selection make tests correlated
and preselected; 20/20 does not imply 20 independent validated mechanisms.

## Limitations
- Small TCGA discovery cohort; residual clinical/technical confounding possible.
- Leave-one-site-out sensitivity analysis was **not performed**.
- Normal samples are adjacent/solid-tissue normals, not healthy donor controls.
- Bulk RNA-seq cannot attribute signatures to specific cell types.
- Cross-platform assay differences and probe mapping may affect replication.
- Candidate genes and pathways were selected using discovery results.
- No survival prediction, clinical utility, or mechanistic causality is established.
- This release audits existing results; it does not rerun all upstream
  GDC downloads and computational analyses in a clean environment.
