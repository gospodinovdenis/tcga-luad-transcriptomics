#!/usr/bin/env python3
"""Release-level consistency checks against shipped results."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
summary = json.loads((TAB/"validation_summary.json").read_text())
genes = pd.read_csv(TAB/"candidate_gene_validation.csv")
pathways = pd.read_csv(TAB/"pathway_signature_validation.csv")
metadata = pd.read_csv(ROOT/"data/metadata/tcga_cohort_metadata.csv")

assert len(genes)==8, "Expected eight preselected candidate genes"
assert genes["external_replication"].astype(str).str.lower().eq("true").all()
assert len(pathways)==20, "Expected twenty preselected pathways"
assert pathways["external_replication"].astype(str).str.lower().eq("true").all()
assert summary["pathways_replicated"] == int(pathways["external_replication"].sum())
assert summary["matched_genes"] == 11431
assert summary["tcga_significant_mapped"] == 2600
assert summary["direction_concordance_n"] == 2401
assert abs(summary["direction_concordance"] - 2401/2600) < 1e-12
assert abs(summary["external_replication_fraction"] - 2061/2600) < 1e-12
assert len(metadata)==30, "Expected 30 discovery samples"
print("PASS: 30 discovery samples, 8 candidates, 20 signatures, validation summary reconciled.")
