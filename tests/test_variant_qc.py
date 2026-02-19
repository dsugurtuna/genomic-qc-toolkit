"""Tests for VariantQC and ConcordanceChecker."""

from genomic_qc.variant_qc import VariantQC, VariantReport
from genomic_qc.concordance import ConcordanceChecker, ConcordanceReport


class TestVariantQC:
    def test_call_rate_check(self):
        qc = VariantQC(min_call_rate=0.98)
        failed = qc.check_call_rates({"rs1": 0.99, "rs2": 0.95, "rs3": 0.98})
        assert failed == ["rs2"]

    def test_batch_effects(self):
        qc = VariantQC(max_batch_diff=0.02)
        failed = qc.check_batch_effects({
            "rs1": {"B1": 0.99, "B2": 0.98},   # diff 0.01 → pass
            "rs2": {"B1": 0.99, "B2": 0.95},   # diff 0.04 → fail
        })
        assert "rs2" in failed
        assert "rs1" not in failed

    def test_find_duplicates(self):
        qc = VariantQC()
        dups = qc.find_duplicates({
            "rs1": "chr1:100",
            "rs2": "chr1:200",
            "rs3": "chr1:100",
        })
        assert "rs1" in dups
        assert "rs3" in dups
        assert "rs2" not in dups

    def test_run_all(self):
        qc = VariantQC()
        report = qc.run(
            call_rates={"rs1": 0.99, "rs2": 0.90},
            batch_rates={},
            positions={"rs1": "chr1:1", "rs2": "chr1:2"},
        )
        assert report.total_variants == 2
        assert len(report.failed_call_rate) == 1


class TestConcordanceChecker:
    def test_compare_sample(self):
        cc = ConcordanceChecker()
        conc, disc, miss = cc.compare_sample(
            {"rs1": "AA", "rs2": "AG", "rs3": "GG"},
            {"rs1": "AA", "rs2": "AG", "rs3": "AA"},
        )
        assert conc == 2
        assert disc == 1
        assert miss == 0

    def test_missing_genotypes(self):
        cc = ConcordanceChecker()
        conc, disc, miss = cc.compare_sample(
            {"rs1": "AA", "rs2": "."},
            {"rs1": "AA", "rs2": "AG"},
        )
        assert miss == 1
        assert conc == 1

    def test_check_flags_low_concordance(self):
        cc = ConcordanceChecker(min_concordance=0.95, min_overlap=2)
        report = cc.check(
            {"S1": {"rs1": "AA", "rs2": "AG", "rs3": "GG"}},
            {"S1": {"rs1": "AA", "rs2": "GG", "rs3": "AA"}},
        )
        assert "S1" in report.flagged_samples
        assert report.concordance_rate < 0.95

    def test_high_concordance(self):
        cc = ConcordanceChecker(min_concordance=0.95, min_overlap=2)
        report = cc.check(
            {"S1": {"rs1": "AA", "rs2": "AG"}},
            {"S1": {"rs1": "AA", "rs2": "AG"}},
        )
        assert report.flagged_samples == []
        assert report.concordance_rate == 1.0
