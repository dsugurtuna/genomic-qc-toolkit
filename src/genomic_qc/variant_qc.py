"""Variant-level QC module.

Flags problematic variants by call rate, missingness per batch,
and genotyping cluster quality.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set


@dataclass
class VariantReport:
    """Variant QC results."""

    total_variants: int = 0
    passed_variants: int = 0
    failed_call_rate: List[str] = field(default_factory=list)
    failed_batch_effect: List[str] = field(default_factory=list)
    failed_duplicate: List[str] = field(default_factory=list)

    @property
    def pass_rate(self) -> float:
        if self.total_variants == 0:
            return 0.0
        return self.passed_variants / self.total_variants

    @property
    def all_failed(self) -> Set[str]:
        return set(self.failed_call_rate) | set(self.failed_batch_effect) | set(self.failed_duplicate)


class VariantQC:
    """Variant-level quality controller.

    Parameters
    ----------
    min_call_rate : float
        Minimum genotyping call rate per variant (default 0.98).
    max_batch_diff : float
        Maximum per-batch missingness difference before flagging (default 0.02).
    """

    def __init__(
        self,
        min_call_rate: float = 0.98,
        max_batch_diff: float = 0.02,
    ) -> None:
        self.min_call_rate = min_call_rate
        self.max_batch_diff = max_batch_diff

    def check_call_rates(
        self, variant_call_rates: Dict[str, float]
    ) -> List[str]:
        """Return variant IDs with call rate below threshold."""
        return [
            vid for vid, cr in variant_call_rates.items()
            if cr < self.min_call_rate
        ]

    def check_batch_effects(
        self, variant_batch_rates: Dict[str, Dict[str, float]]
    ) -> List[str]:
        """Detect variants with inconsistent call rates across batches.

        Parameters
        ----------
        variant_batch_rates : dict
            {variant_id: {batch_id: call_rate}}
        """
        flagged: List[str] = []
        for vid, batch_rates in variant_batch_rates.items():
            rates = list(batch_rates.values())
            if len(rates) < 2:
                continue
            if max(rates) - min(rates) > self.max_batch_diff:
                flagged.append(vid)
        return flagged

    def find_duplicates(
        self, variant_positions: Dict[str, str]
    ) -> List[str]:
        """Find duplicate variants by chromosomal position.

        Parameters
        ----------
        variant_positions : dict
            {variant_id: "chr:pos"}
        """
        seen: Dict[str, str] = {}
        duplicates: List[str] = []
        for vid, pos in variant_positions.items():
            if pos in seen:
                duplicates.append(vid)
                if seen[pos] not in duplicates:
                    duplicates.append(seen[pos])
            else:
                seen[pos] = vid
        return duplicates

    def run(
        self,
        call_rates: Dict[str, float],
        batch_rates: Dict[str, Dict[str, float]],
        positions: Dict[str, str],
    ) -> VariantReport:
        """Run all variant QC checks."""
        report = VariantReport(total_variants=len(call_rates))
        report.failed_call_rate = self.check_call_rates(call_rates)
        report.failed_batch_effect = self.check_batch_effects(batch_rates)
        report.failed_duplicate = self.find_duplicates(positions)
        report.passed_variants = report.total_variants - len(report.all_failed)
        return report
