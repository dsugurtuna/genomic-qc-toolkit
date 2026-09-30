"""Genotype concordance between two platforms (for example array and WGS).

Low concordance for one sample usually means a sample swap or
contamination; low concordance for everyone usually means a strand,
build or reference-allele mismatch between the two call sets.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

MISSING = frozenset({"", ".", "./.", ".|.", "NA", "N/A", "--", "0 0", "00", "NN"})


def normalise_genotype(genotype: str) -> str | None:
    """Return an unordered genotype such as ``A/G``, or None if missing.

    ``AG``, ``G/A``, ``A|G`` and ``A G`` all become ``A/G``. VCF GT codes are
    kept as codes (``0|1`` becomes ``0/1``); note that ``0/0`` is homozygous
    reference, not missing.
    """
    text = genotype.strip().upper()
    if text in MISSING:
        return None
    alleles = re.split(r"[/| ]", text) if re.search(r"[/| ]", text) else list(text)
    if len(alleles) != 2 or any(a in ("", ".") for a in alleles):
        return None
    return "/".join(sorted(alleles))


@dataclass
class ConcordanceReport:
    """Concordance across all samples present on both platforms."""

    total_comparisons: int = 0  # samples compared
    concordant: int = 0
    discordant: int = 0
    missing: int = 0
    flagged_samples: list[str] = field(default_factory=list)
    insufficient_overlap: list[str] = field(default_factory=list)
    per_sample_rate: dict[str, float] = field(default_factory=dict)

    @property
    def concordance_rate(self) -> float:
        total = self.concordant + self.discordant
        return self.concordant / total if total else 0.0


class ConcordanceChecker:
    """Compare genotype calls between two platforms.

    Parameters
    ----------
    min_concordance : float
        Per-sample concordance below which a sample is flagged (default 0.98).
    min_overlap : int
        Non-missing shared variants needed to judge a sample (default 100).
        Samples with fewer are listed in ``insufficient_overlap``.
    """

    def __init__(self, min_concordance: float = 0.98, min_overlap: int = 100) -> None:
        self.min_concordance = min_concordance
        self.min_overlap = min_overlap

    def compare_sample(
        self, calls_a: dict[str, str], calls_b: dict[str, str]
    ) -> tuple[int, int, int]:
        """Return (concordant, discordant, missing) over shared variants."""
        concordant = discordant = missing = 0
        for vid in calls_a.keys() & calls_b.keys():
            ga, gb = normalise_genotype(calls_a[vid]), normalise_genotype(calls_b[vid])
            if ga is None or gb is None:
                missing += 1
            elif ga == gb:
                concordant += 1
            else:
                discordant += 1
        return concordant, discordant, missing

    def check(
        self,
        platform_a: dict[str, dict[str, str]],
        platform_b: dict[str, dict[str, str]],
    ) -> ConcordanceReport:
        """Check every sample present on both platforms.

        Parameters
        ----------
        platform_a, platform_b : dict
            {sample_id: {variant_id: genotype}}
        """
        report = ConcordanceReport()
        shared = sorted(platform_a.keys() & platform_b.keys())
        report.total_comparisons = len(shared)
        for sid in shared:
            conc, disc, miss = self.compare_sample(platform_a[sid], platform_b[sid])
            report.concordant += conc
            report.discordant += disc
            report.missing += miss
            if conc + disc < self.min_overlap:
                report.insufficient_overlap.append(sid)
                continue
            rate = conc / (conc + disc)
            report.per_sample_rate[sid] = rate
            if rate < self.min_concordance:
                report.flagged_samples.append(sid)
        return report
