import subprocess, sys
from pathlib import Path

def test_release_audit():
    root=Path(__file__).resolve().parents[1]
    result=subprocess.run([sys.executable,str(root/"scripts/audit_release.py")],
                          cwd=root,capture_output=True,text=True)
    assert result.returncode==0, result.stdout+"\n"+result.stderr

def test_required_outputs():
    root=Path(__file__).resolve().parents[1]
    for name in ["README.md","docs/METHODS_AND_LIMITATIONS.md",
                 "results/tables/validation_summary.json",
                 "results/figures/cross_cohort_concordance.png"]:
        assert (root/name).is_file()
