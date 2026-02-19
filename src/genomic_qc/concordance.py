"""Concordance checker module.

Cross-validates genotype calls between platforms (e.g. array vs WGS)
to detect sample swaps or systematic calling errors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple


@dataclass
class ConcordanceReport:
    """Concordance check results."""

    total_comparisons: int = 0
    concordant: int = 0
    discordant: int = 0
    missing: int = 0
    flagged_samples: List[str] = field(default_factory=list)

    @property
    def concordance_rate(self) -> float:
        total = self.concordant + self.discordant
        if total == 0:
            return 0.0
        return self.concordant / total


class ConcordanceChecker:
    """Compare genotype calls between two platforms.

    Parameters
    ----------
    min_concordance : float
        Minimum per-sample concordance rate (default 0.98).
    min_overlap : int
        Minimum overlapping variants required (default 100).
    """

    def __init__(
        self,
        min_concordance: float = 0.98,
        min_overlap: int = 100,
    ) -> None:
        self.min_concordance = min_concordance
        self.min_overlap = min_overlap

    def compare_sample(
        self,
        calls_a: Dict[str, str],
        calls_b: Dict[str, str],
    ) -> Tuple[int, int, int]:
        """Compare genotype calls for a single sample.

        Parameters
        ----------
        calls_a, calls_b : dict
            {variant_id: genotype_string}

        Returns
        -------
        (concordant, discordant, missing)
        """
        overlap = set(calls_a.keys()) & set(calls_b.keys())
        concordant = 0
        discordant = 0
        missing = 0
        for vid in overlap:
            ga, gb = calls_a[vid], calls_b[vid]
            if ga == "." or gb == ".":
                missing += 1
            elif ga == gb:
                concordant += 1
            else:
                discordant += 1
        return concordant, discordant, missing

    def check(
        self,
        platform_a: Dict[str, Dict[str, str]],
        platform_b: Dict[str, Dict[str, str]],
    ) -> ConcordanceReport:
        """Check concordance across all shared samples.

        Parameters
        ----------
        platform_a, platform_b : dict
            {sample_id: {variant_id: genotype}}
        """
        report = ConcordanceReport()
        shared = set(platform_a.keys()) & set(platform_b.keys())
        report.total_comparisons = len(shared)

        for sid in shared:
            conc, disc, miss = self.compare_sample(
                platform_a[sid], platform_b[sid]
            )
            report.concordant += conc
            report.discordant += disc
            report.missing += miss
            total = conc + disc
            if total >= self.min_overlap:
                rate = conc / total
                if rate < self.min_concordance:
                    report.flagged_samples.append(sid)
        return report
