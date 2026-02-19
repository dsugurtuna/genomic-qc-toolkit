"""Sample-level QC module.

Identifies problematic samples via contamination estimation,
sex verification, and coverage checks.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set


@dataclass
class SampleReport:
    """Sample QC results."""

    total_samples: int = 0
    passed_samples: int = 0
    failed_contamination: List[str] = field(default_factory=list)
    failed_sex_check: List[str] = field(default_factory=list)
    failed_coverage: List[str] = field(default_factory=list)

    @property
    def pass_rate(self) -> float:
        if self.total_samples == 0:
            return 0.0
        return self.passed_samples / self.total_samples

    @property
    def all_failed(self) -> Set[str]:
        return set(self.failed_contamination) | set(self.failed_sex_check) | set(self.failed_coverage)


class SampleQC:
    """Sample-level quality controller.

    Parameters
    ----------
    max_contamination : float
        Maximum acceptable contamination fraction (default 0.05).
    min_coverage : float
        Minimum acceptable mean coverage (default 20.0).
    f_threshold_male : float
        Minimum X-chromosome F statistic to classify as male (default 0.8).
    f_threshold_female : float
        Maximum X-chromosome F statistic to classify as female (default 0.2).
    """

    def __init__(
        self,
        max_contamination: float = 0.05,
        min_coverage: float = 20.0,
        f_threshold_male: float = 0.8,
        f_threshold_female: float = 0.2,
    ) -> None:
        self.max_contamination = max_contamination
        self.min_coverage = min_coverage
        self.f_threshold_male = f_threshold_male
        self.f_threshold_female = f_threshold_female

    def check_contamination(
        self, sample_contamination: Dict[str, float]
    ) -> List[str]:
        """Return sample IDs with contamination above threshold."""
        return [
            sid for sid, cont in sample_contamination.items()
            if cont > self.max_contamination
        ]

    def check_sex(
        self, sample_sex: Dict[str, tuple]
    ) -> List[str]:
        """Check sex assignment concordance.

        Parameters
        ----------
        sample_sex : dict
            {sample_id: (reported_sex, f_statistic)} where
            reported_sex is 1 (male) or 2 (female).
        """
        failed: List[str] = []
        for sid, (reported, f_stat) in sample_sex.items():
            if reported == 1 and f_stat < self.f_threshold_male:
                failed.append(sid)
            elif reported == 2 and f_stat > self.f_threshold_female:
                failed.append(sid)
        return failed

    def check_coverage(
        self, sample_coverage: Dict[str, float]
    ) -> List[str]:
        """Return samples below minimum coverage."""
        return [
            sid for sid, cov in sample_coverage.items()
            if cov < self.min_coverage
        ]

    def run(
        self,
        contamination: Dict[str, float],
        sex: Dict[str, tuple],
        coverage: Dict[str, float],
    ) -> SampleReport:
        """Run all sample QC checks."""
        all_samples = set(contamination.keys()) | set(sex.keys()) | set(coverage.keys())
        report = SampleReport(total_samples=len(all_samples))
        report.failed_contamination = self.check_contamination(contamination)
        report.failed_sex_check = self.check_sex(sex)
        report.failed_coverage = self.check_coverage(coverage)
        report.passed_samples = report.total_samples - len(report.all_failed)
        return report
