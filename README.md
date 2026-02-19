# Genomic QC Toolkit

[![CI](https://github.com/dsugurtuna/genomic-qc-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/dsugurtuna/genomic-qc-toolkit/actions/workflows/ci.yml)

**Sample and variant quality control for WGS, WES, and genotyping array data.**

Provides a structured QC pipeline covering sample-level checks (contamination, sex verification, coverage), variant-level checks (call rate, batch effects, duplicate positions), and cross-platform concordance validation.

> **Portfolio project.** Demonstrates generalised genomic QC workflows. No real participant data is included.

---

## Architecture

```
src/genomic_qc/
    __init__.py        # Public API exports
    sample_qc.py       # Sample QC (contamination, sex, coverage)
    variant_qc.py      # Variant QC (call rate, batch effects, duplicates)
    concordance.py     # Cross-platform genotype concordance
tests/
    test_sample_qc.py  # Sample QC tests
    test_variant_qc.py # Variant QC and concordance tests
```

---

## Quick start

```bash
pip install -e ".[dev]"
pytest -v
```

### Python API

```python
from genomic_qc import SampleQC, VariantQC, ConcordanceChecker

# Sample QC
sqc = SampleQC(max_contamination=0.05, min_coverage=20.0)
report = sqc.run(
    contamination={"S001": 0.01, "S002": 0.12},
    sex={"S001": (1, 0.95), "S002": (2, 0.10)},
    coverage={"S001": 30.5, "S002": 18.0},
)
print(f"Sample pass rate: {report.pass_rate:.1%}")

# Variant QC
vqc = VariantQC(min_call_rate=0.98, max_batch_diff=0.02)
vreport = vqc.run(
    call_rates={"rs1": 0.99, "rs2": 0.90},
    batch_rates={"rs1": {"B1": 0.99, "B2": 0.98}},
    positions={"rs1": "chr1:100", "rs2": "chr1:200"},
)

# Cross-platform concordance
cc = ConcordanceChecker(min_concordance=0.98)
conc_report = cc.check(array_data, wgs_data)
print(f"Concordance: {conc_report.concordance_rate:.1%}")
```

---

## Key features

| Feature | Detail |
| :--- | :--- |
| **Contamination check** | Flags samples above configurable FreeMix threshold |
| **Sex verification** | X-chromosome F-statistic concordance with reported sex |
| **Coverage check** | Mean depth filtering for WGS/WES samples |
| **Variant call rate** | Per-variant missingness filtering |
| **Batch effect detection** | Flags variants with inconsistent missingness across batches |
| **Duplicate detection** | Identifies variants sharing chromosomal position |
| **Cross-platform concordance** | Genotype comparison between array and sequencing platforms |

## Development

```bash
make dev        # install with dev dependencies
make test       # run pytest
make lint       # run ruff
make clean      # remove build artefacts
```

## Jira provenance

| Ticket | Description |
| :--- | :--- |
| BIOIN-684 | End-to-end QC pipeline for genotyping array and WGS/WES data |

---

*Created by [dsugurtuna](https://github.com/dsugurtuna)*
