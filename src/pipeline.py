"""Public GDC TCGA-LUAD STAR counts pilot workflow (no authentication)."""
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://api.gdc.cancer.gov"
MEASURES = {"unstranded": "counts", "tpm_unstranded": "tpm", "fpkm_unstranded": "fpkm"}
META_ROWS = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}


def session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET", "POST"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": "tcga-luad-colab-workflow/1.0"})
    return s


def fetch_manifest(s, page_size=500):
    filters = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ]}
    fields = "file_id,file_name,md5sum,file_size,data_type,analysis.workflow_type,cases.case_id,cases.submitter_id,cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type"
    records = []
    offset = 0
    while True:
        params = {"filters": json.dumps(filters), "fields": fields, "format": "JSON", "size": page_size, "from": offset}
        resp = s.get(f"{BASE}/files", params=params, timeout=120)
        resp.raise_for_status()
        data = resp.json()["data"]
        hits = data["hits"]
        records.extend(hits)
        offset += len(hits)
        if not hits or offset >= int(data["pagination"]["total"]):
            break
    if not records:
        raise RuntimeError("No matching open-access TCGA-LUAD STAR counts files were found")
    return records


def flatten_manifest(records):
    out = []
    for f in records:
        cases = f.get("cases") or []
        # Do not silently assign a file to an ambiguous case/sample.
        if len(cases) != 1 or len(cases[0].get("samples") or []) != 1:
            continue
        c = cases[0]
        sample = c["samples"][0]
        out.append({"file_id": f["file_id"], "file_name": f.get("file_name"), "md5sum": f.get("md5sum"),
                    "file_size": f.get("file_size"), "case_id": c.get("case_id"),
                    "case_barcode": c.get("submitter_id"), "sample_id": sample.get("sample_id"),
                    "sample_barcode": sample.get("submitter_id"), "sample_type": sample.get("sample_type")})
    df = pd.DataFrame(out)
    if df.empty:
        raise RuntimeError("No unambiguous file-to-sample records; inspect the GDC manifest")
    if df.file_id.duplicated().any():
        raise ValueError("Duplicate file UUIDs in manifest")
    return df.sort_values(["sample_type", "sample_barcode", "file_id"]).reset_index(drop=True)


def select_pilot(manifest, tumor_n=8, normal_n=4):
    chosen = []
    for kind, n in [("Primary Tumor", tumor_n), ("Solid Tissue Normal", normal_n)]:
        subset = manifest[manifest.sample_type.eq(kind)].sort_values(["sample_barcode", "file_id"])
        # Avoid replicate files for the same sample in the pilot.
        subset = subset.drop_duplicates("sample_barcode")
        chosen.append(subset.head(n))
        if len(subset) < n:
            print(f"Warning: requested {n} {kind} samples, found {len(subset)}")
    pilot = pd.concat(chosen, ignore_index=True)
    if pilot.empty:
        raise RuntimeError("No tumor/normal samples available")
    if pilot.sample_barcode.isna().any() or pilot.sample_barcode.duplicated().any():
        raise ValueError("Missing or duplicate sample barcodes in selected pilot")
    return pilot


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download_file(s, file_id, destination, expected_md5=None):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and (not expected_md5 or md5_file(destination).lower() == expected_md5.lower()):
        return destination
    tmp = destination.with_suffix(destination.suffix + ".part")
    try:
        with s.get(f"{BASE}/data/{file_id}", stream=True, timeout=(30, 180)) as r:
            r.raise_for_status()
            with open(tmp, "wb") as out:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        out.write(chunk)
        if expected_md5 and md5_file(tmp).lower() != expected_md5.lower():
            raise ValueError(f"MD5 mismatch: {file_id}")
        tmp.replace(destination)
    finally:
        tmp.unlink(missing_ok=True)
    return destination


def parse_star_counts(path):
    # STAR counts files have a leading comment line and 4 STAR summary rows.
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"gene_id": "string"}, low_memory=False)
    required = {"gene_id", "gene_name", "gene_type", *MEASURES}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing expected STAR columns in {path}: {sorted(missing)}; found {list(df.columns)}")
    df = df.loc[~df["gene_id"].isin(META_ROWS)].copy()
    if df.empty or df.gene_id.isna().any() or df.gene_id.duplicated().any():
        raise ValueError(f"Invalid gene IDs in {path}")
    for measure in MEASURES:
        df[measure] = pd.to_numeric(df[measure], errors="raise")
        if not np.isfinite(df[measure].to_numpy(dtype=float)).all() or (df[measure] < 0).any():
            raise ValueError(f"Non-finite or negative {measure} values in {path}")
    if not (df["unstranded"] == np.floor(df["unstranded"])).all():
        raise ValueError(f"Non-integer unstranded counts in {path}")
    return df


def assemble(pilot, raw_dir):
    count_series, tpm_series, fpkm_series = [], [], []
    annotations = None
    for row in pilot.itertuples(index=False):
        frame = parse_star_counts(Path(raw_dir) / f"{row.file_id}.tsv")
        ann = frame[["gene_id", "gene_name", "gene_type"]].set_index("gene_id")
        if annotations is None:
            annotations = ann
        elif not annotations.index.equals(ann.index):
            raise ValueError(f"Gene identifiers/order mismatch in {row.file_id}; resolve before merging")
        elif not annotations.equals(ann):
            raise ValueError(f"Gene annotation mismatch in {row.file_id}; resolve before merging")
        for field, target in [("unstranded", count_series), ("tpm_unstranded", tpm_series), ("fpkm_unstranded", fpkm_series)]:
            target.append(frame.set_index("gene_id")[field].rename(row.sample_barcode))
    matrices = {"counts": pd.concat(count_series, axis=1), "tpm": pd.concat(tpm_series, axis=1), "fpkm": pd.concat(fpkm_series, axis=1)}
    validate_matrices(matrices, pilot)
    return matrices, annotations


def validate_matrices(matrices, metadata):
    counts = matrices["counts"]
    if counts.empty or counts.index.has_duplicates or counts.columns.has_duplicates:
        raise ValueError("Empty/duplicate count matrix")
    if list(counts.columns) != metadata.sample_barcode.tolist():
        raise ValueError("Metadata/expression sample ordering mismatch")
    for label, matrix in matrices.items():
        if not matrix.index.equals(counts.index) or not matrix.columns.equals(counts.columns):
            raise ValueError(f"Matrix axes differ: {label}")
        arr = matrix.to_numpy(dtype=float)
        if not np.isfinite(arr).all() or (arr < 0).any():
            raise ValueError(f"Invalid values in {label}")
        if label == "counts" and not np.equal(arr, np.floor(arr)).all():
            raise ValueError("Raw counts must be integers")
    return True


def qc_report(matrices, annotations, metadata, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts = matrices["counts"]
    totals = counts.sum(axis=0)
    detected = (counts > 0).sum(axis=0)
    sample_qc = metadata.set_index("sample_barcode").loc[counts.columns, ["case_barcode", "sample_type"]].copy()
    sample_qc["total_assigned_gene_counts"] = totals
    sample_qc["detected_genes"] = detected
    sample_qc["zero_fraction"] = (counts == 0).mean(axis=0)
    sample_qc["low_count_total_flag"] = totals < totals.median() * 0.25
    sample_qc.to_csv(out / "sample_qc.csv")
    fig, ax = plt.subplots(figsize=(max(8, len(counts.columns) * 0.65), 4))
    ax.bar(range(len(totals)), totals.to_numpy() / 1e6)
    ax.set_xticks(range(len(totals)), counts.columns, rotation=90, fontsize=7)
    ax.set_ylabel("Assigned gene counts (millions)")
    ax.set_title("TCGA-LUAD STAR unstranded counts: sample totals")
    fig.tight_layout()
    fig.savefig(out / "library_size.png", dpi=160)
    plt.close(fig)
    summary = {"n_samples": int(counts.shape[1]), "n_genes": int(counts.shape[0]),
               "n_primary_tumor": int(metadata.sample_type.eq("Primary Tumor").sum()),
               "n_solid_tissue_normal": int(metadata.sample_type.eq("Solid Tissue Normal").sum()),
               "n_flagged_low_count": int(sample_qc.low_count_total_flag.sum()),
               "missing_gene_names": int(annotations.gene_name.isna().sum()),
               "validation": "PASS", "note": "Expression-level checks only; not FASTQ/read-level sequencing QC."}
    (out / "qc_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def run(output_dir="tcga_luad_output", tumor_n=8, normal_n=4):
    out = Path(output_dir)
    for sub in ["raw", "processed", "qc", "provenance"]:
        (out / sub).mkdir(parents=True, exist_ok=True)
    s = session()
    records = fetch_manifest(s)
    (out / "provenance" / "gdc_manifest_full.json").write_text(json.dumps(records, indent=2))
    manifest = flatten_manifest(records)
    manifest.to_csv(out / "provenance" / "gdc_manifest_flat.csv", index=False)
    pilot = select_pilot(manifest, tumor_n, normal_n)
    pilot.to_csv(out / "processed" / "sample_metadata.csv", index=False)
    for i, row in enumerate(pilot.itertuples(index=False), 1):
        print(f"Downloading {i}/{len(pilot)}: {row.sample_barcode} ({row.file_id})", flush=True)
        download_file(s, row.file_id, out / "raw" / f"{row.file_id}.tsv", row.md5sum)
    matrices, annotations = assemble(pilot, out / "raw")
    for name, matrix in matrices.items():
        matrix.to_csv(out / "processed" / f"{name}_matrix.csv.gz", compression="gzip")
    annotations.to_csv(out / "processed" / "gene_annotations.csv")
    summary = qc_report(matrices, annotations, pilot, out / "qc")
    versions = {"python": sys.version, "platform": platform.platform(), "pandas": pd.__version__,
                "numpy": np.__version__, "requests": requests.__version__, "matplotlib": matplotlib.__version__,
                "gdc_api": BASE, "cohort": "TCGA-LUAD", "workflow": "STAR - Counts", "selection": {"tumor_n": tumor_n, "normal_n": normal_n}}
    (out / "provenance" / "run_environment.json").write_text(json.dumps(versions, indent=2))
    print("QC summary:", json.dumps(summary, indent=2))
    return summary
