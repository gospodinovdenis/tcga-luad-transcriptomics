# Reproduction guide

## Quick audit (does not require downloading sequencing data)
```bash
python -m pip install -r requirements.txt
python scripts/audit_release.py
python -m pytest -q tests/test_release.py
```

## Full upstream analyses
The notebooks in `notebooks/` are preserved from the original phase
archives and retain their original **archive upload/path conventions**.
They are **historical execution notebooks**, not a unified one-click
pipeline. To rerun, use the relevant original phase ZIP and follow the
instructions in that phase's notebook. Some require Google Colab, R,
Bioconductor, GEO/GDC internet access, and external input archives.

1. Phase 1: GDC download and QC.
2. Phase 4.5: pilot QC.
3. Phase 5: pilot DESeq2 normalization.
4. Phase 6: site-balanced discovery cohort and DESeq2.
5. Phase 7: volcano and heatmap.
6. Phase 7B: GO/Reactome ORA and preranked GSEA.
7. Phase 8: interpretation report.
8. Phase 9: independent GEO paired validation.

The curated result tables here permit an independent **result-level audit**
without rerunning the full workflow. Do not claim end-to-end clean-room
reproduction until a complete rerun has been tested.
