# Data sources and provenance

## Discovery: TCGA-LUAD
National Cancer Institute Genomic Data Commons (GDC), TCGA-LUAD project,
STAR - Counts workflow, raw unstranded integer counts. Discovery cohort:
30 samples (15 primary tumor, 15 solid tissue normal), seven overlapping
tissue-source sites. Sample IDs are in `data/metadata/tcga_cohort_metadata.csv`.
The site-adjusted DESeq2 contrast is tumor vs normal, with `~ tss + sample_type`.
The original pilot (12 samples) was replaced by the site-balanced exploratory
cohort for inference.

GDC portal: https://portal.gdc.cancer.gov/projects/TCGA-LUAD

## External validation: GEO GSE32863
NCBI GEO accession GSE32863, platform GPL6884, Illumina HumanWG-6 v3.0
microarray. 116 samples total: 58 tumor and 58 adjacent non-tumor lung.
57 complete tumor/normal pairs retained. Patient 3023 has only a tumor
sample; patient 3035 has only a normal sample. These two samples are excluded
from paired validation.

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32863

The external study is **not** healthy-donor lung tissue. GEO data are
retrieved on demand in the Phase 9 notebook. Source-provided normalized
log2 microarray intensities were analyzed using limma with patient fixed
effects. Probe selection used maximum across-sample variance per gene
symbol, independent of phenotype. Cross-platform logFC magnitudes are
not directly interchangeable.

## Redistribution
The repository distributes small derived result tables and public study
metadata. It does not redistribute raw GDC downloads or full external
microarray expression matrices. Users must comply with upstream data
source terms and cite the original studies.
