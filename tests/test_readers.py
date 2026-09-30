"""Tests for the tool-output readers and the sample QC demo."""

import runpy
from pathlib import Path

import pytest

from genomic_qc.readers import (
    read_lmiss,
    read_mosdepth_summaries,
    read_selfsm,
    read_sexcheck,
)

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_readers_on_examples() -> None:
    assert read_selfsm(EXAMPLES / "cohort.selfSM")["SYN03"] == 0.083
    assert read_sexcheck(EXAMPLES / "cohort.sexcheck")["SYN06"] == (2, 0.97)
    depth = read_mosdepth_summaries(sorted((EXAMPLES / "mosdepth").glob("*.txt")))
    assert depth["SYN08"] == 14.2
    assert len(depth) == 8


def test_read_lmiss(tmp_path: Path) -> None:
    lmiss = tmp_path / "x.lmiss"
    lmiss.write_text(" CHR  SNP  N_MISS  N_GENO  F_MISS\n   1  rs1  1  100  0.01\n")
    assert read_lmiss(lmiss) == {"rs1": pytest.approx(0.99)}


def test_sample_demo(capsys: pytest.CaptureFixture[str]) -> None:
    runpy.run_path(str(EXAMPLES / "sample_qc_demo.py"), run_name="__main__")
    out = capsys.readouterr().out
    assert "contamination > 5%: ['SYN03']" in out
    assert "sex mismatch:       ['SYN06']" in out
    assert "coverage < 20x:     ['SYN08']" in out
    assert "passed 5 of 8" in out
