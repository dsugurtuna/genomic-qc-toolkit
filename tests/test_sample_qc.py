"""Tests for SampleQC."""

from genomic_qc.sample_qc import SampleQC, SampleReport


class TestSampleQC:
    def test_contamination_check(self):
        qc = SampleQC(max_contamination=0.05)
        failed = qc.check_contamination({"S1": 0.01, "S2": 0.10, "S3": 0.04})
        assert failed == ["S2"]

    def test_sex_check(self):
        qc = SampleQC()
        failed = qc.check_sex({
            "S1": (1, 0.95),   # male, high F → pass
            "S2": (1, 0.30),   # male, low F → fail
            "S3": (2, 0.10),   # female, low F → pass
            "S4": (2, 0.85),   # female, high F → fail
        })
        assert "S2" in failed
        assert "S4" in failed
        assert "S1" not in failed
        assert "S3" not in failed

    def test_coverage_check(self):
        qc = SampleQC(min_coverage=20.0)
        failed = qc.check_coverage({"S1": 30.0, "S2": 15.0, "S3": 20.0})
        assert failed == ["S2"]

    def test_run_all(self):
        qc = SampleQC(max_contamination=0.05, min_coverage=20.0)
        report = qc.run(
            contamination={"S1": 0.01, "S2": 0.10},
            sex={"S1": (1, 0.95), "S2": (1, 0.95)},
            coverage={"S1": 30.0, "S2": 25.0},
        )
        assert report.total_samples == 2
        assert report.passed_samples == 1
        assert len(report.failed_contamination) == 1

    def test_pass_rate(self):
        r = SampleReport(total_samples=100, passed_samples=95)
        assert r.pass_rate == 0.95

    def test_empty_pass_rate(self):
        r = SampleReport()
        assert r.pass_rate == 0.0
