"""Variant-level QC thresholds.

Flags variants by overall call rate, by the spread of call rates between
batches (a simple batch-effect screen) and by duplicate positions.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class VariantReport:
    """Variant QC results."""

    total_variants: int = 0
    passed_variants: int = 0
    failed_call_rate: list[str] = field(default_factory=list)
    failed_batch_effect: list[str] = field(default_factory=list)
    failed_duplicate: list[str] = field(default_factory=list)

    @property
    def pass_rate(self) -> float:
        if self.total_variants == 0:
            return 0.0
        return self.passed_variants / self.total_variants

    @property
    def all_failed(self) -> set[str]:
        return (
            set(self.failed_call_rate)
            | set(self.failed_batch_effect)
            | set(self.failed_duplicate)
        )


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

    def check_call_rates(self, variant_call_rates: dict[str, float]) -> list[str]:
        """Return variant IDs with call rate below threshold."""
        return [
            vid for vid, cr in variant_call_rates.items() if cr < self.min_call_rate
        ]

    def check_batch_effects(
        self, variant_batch_rates: dict[str, dict[str, float]]
    ) -> list[str]:
        """Detect variants with inconsistent call rates across batches.

        Parameters
        ----------
        variant_batch_rates : dict
            {variant_id: {batch_id: call_rate}}
        """
        flagged: list[str] = []
        for vid, batch_rates in variant_batch_rates.items():
            rates = list(batch_rates.values())
            if len(rates) < 2:
                continue
            # A range, not a statistical test: a quick screen for variants
            # that behave differently in one batch.
            if max(rates) - min(rates) > self.max_batch_diff:
                flagged.append(vid)
        return flagged

    @staticmethod
    def _normalise_position(pos: str) -> str:
        """``chr1:100`` and ``1:100`` are the same position."""
        pos = pos.strip()
        return pos[3:] if pos.lower().startswith("chr") else pos

    def find_duplicates(self, variant_positions: dict[str, str]) -> list[str]:
        """Variants that share a position with another variant.

        Parameters
        ----------
        variant_positions : dict
            {variant_id: "chr:pos"}; a ``chr`` prefix is ignored.
        """
        by_position: dict[str, list[str]] = {}
        for vid, pos in variant_positions.items():
            by_position.setdefault(self._normalise_position(pos), []).append(vid)
        return [vid for ids in by_position.values() if len(ids) > 1 for vid in ids]

    def run(
        self,
        call_rates: dict[str, float],
        batch_rates: dict[str, dict[str, float]],
        positions: dict[str, str],
    ) -> VariantReport:
        """Run all variant QC checks.

        The total counts every variant seen in any input, so a variant that
        appears only in ``positions`` still counts towards the pass rate.
        """
        all_variants = set(call_rates) | set(batch_rates) | set(positions)
        report = VariantReport(total_variants=len(all_variants))
        report.failed_call_rate = self.check_call_rates(call_rates)
        report.failed_batch_effect = self.check_batch_effects(batch_rates)
        report.failed_duplicate = self.find_duplicates(positions)
        report.passed_variants = report.total_variants - len(report.all_failed)
        return report
