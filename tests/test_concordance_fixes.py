"""Tests for concordance normalisation, duplicate positions and variant totals."""

import pytest

from genomic_qc.concordance import ConcordanceChecker, normalise_genotype
from genomic_qc.variant_qc import VariantQC


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("AG", "A/G"),
        ("G/A", "A/G"),
        ("g|a", "A/G"),
        ("A G", "A/G"),
        ("0|1", "0/1"),
        ("0/0", "0/0"),  # homozygous reference, not missing
        ("./.", None),
        ("--", None),
        ("0 0", None),  # PLINK missing
        ("NA", None),
    ],
)
def test_normalise_genotype(raw: str, expected: str | None) -> None:
    assert normalise_genotype(raw) == expected


def test_concordance_ignores_allele_order() -> None:
    cc = ConcordanceChecker(min_overlap=1)
    assert cc.compare_sample(
        {"rs1": "AG", "rs2": "CC"}, {"rs1": "G/A", "rs2": "CC"}
    ) == (
        2,
        0,
        0,
    )


def test_insufficient_overlap_is_reported_not_ignored() -> None:
    cc = ConcordanceChecker(min_concordance=0.9, min_overlap=3)
    report = cc.check(
        {"S1": {"a": "AA", "b": "AG", "c": "GG"}, "S2": {"a": "AA"}},
        {"S1": {"a": "AA", "b": "AA", "c": "AA"}, "S2": {"a": "AA"}},
    )
    assert report.flagged_samples == ["S1"]
    assert report.per_sample_rate["S1"] == pytest.approx(1 / 3)
    assert report.insufficient_overlap == ["S2"]


def test_duplicates_ignore_chr_prefix() -> None:
    dups = VariantQC().find_duplicates(
        {"rs1": "chr1:100", "rs2": "1:100", "rs3": "1:200"}
    )
    assert sorted(dups) == ["rs1", "rs2"]


def test_variant_total_counts_every_input() -> None:
    report = VariantQC().run(
        call_rates={"rs1": 0.99},
        batch_rates={},
        positions={"rs1": "1:100", "rs2": "1:100"},
    )
    assert report.total_variants == 2
    assert report.passed_variants == 0
