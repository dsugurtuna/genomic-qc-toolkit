"""Sample QC on the synthetic tool outputs in this folder.

python examples/sample_qc_demo.py
"""

from pathlib import Path

from genomic_qc import SampleQC
from genomic_qc.readers import read_mosdepth_summaries, read_selfsm, read_sexcheck

HERE = Path(__file__).parent

report = SampleQC(max_contamination=0.05, min_coverage=20.0).run(
    contamination=read_selfsm(HERE / "cohort.selfSM"),
    sex=read_sexcheck(HERE / "cohort.sexcheck"),
    coverage=read_mosdepth_summaries(sorted((HERE / "mosdepth").glob("*.txt"))),
)
print(f"contamination > 5%: {report.failed_contamination}")
print(f"sex mismatch:       {report.failed_sex_check}")
print(f"coverage < 20x:     {report.failed_coverage}")
print(f"passed {report.passed_samples} of {report.total_samples}")
