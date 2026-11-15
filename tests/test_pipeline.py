import pandas as pd
import pytest
from src.pipeline import parse_star_counts, validate_matrices, select_pilot


def test_parser(tmp_path):
    f = tmp_path / "star.tsv"
    f.write_text("# gene-model: test\ngene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\nN_unmapped\tNA\tNA\t3\t3\t3\t0\t0\t0\nENSG000001\tA\tprotein_coding\t10\t5\t5\t2.5\t1.2\t1.4\n")
    df = parse_star_counts(f)
    assert len(df) == 1
    assert df.iloc[0].unstranded == 10


def test_invalid_count(tmp_path):
    f = tmp_path / "bad.tsv"
    f.write_text("gene_id\tgene_name\tgene_type\tunstranded\ttpm_unstranded\tfpkm_unstranded\nENSG1\tA\tprotein_coding\t-1\t0\t0\n")
    with pytest.raises(ValueError):
        parse_star_counts(f)


def test_sample_order_validation():
    counts = pd.DataFrame({"S1": [1, 2], "S2": [0, 4]}, index=["g1", "g2"])
    meta = pd.DataFrame({"sample_barcode": ["S2", "S1"]})
    with pytest.raises(ValueError, match="ordering"):
        validate_matrices({"counts": counts, "tpm": counts.astype(float), "fpkm": counts.astype(float)}, meta)


def test_pilot_balanced():
    df = pd.DataFrame({"sample_type": ["Primary Tumor", "Primary Tumor", "Solid Tissue Normal"], "sample_barcode": ["A", "B", "C"], "file_id": ["1", "2", "3"]})
    pilot = select_pilot(df, tumor_n=1, normal_n=1)
    assert len(pilot) == 2
