"""Loaders for common QC tool outputs, returning the dicts the checks expect.

* VerifyBamID / VerifyBamID2 ``.selfSM``: ``FREEMIX`` per ``#SEQ_ID``.
* PLINK ``--check-sex`` ``.sexcheck``: ``(PEDSEX, F)`` per ``IID``.
* mosdepth ``<sample>.mosdepth.summary.txt``: the ``mean`` of the
  ``total`` row (whole-genome mean depth).
* PLINK ``--missing`` ``.lmiss``: call rate (1 - F_MISS) per ``SNP``.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path


def _table(path: str | Path) -> Iterator[dict[str, str]]:
    with open(path) as fh:
        header = fh.readline().split()
        for line in fh:
            parts = line.split()
            if len(parts) == len(header):
                yield dict(zip(header, parts, strict=True))


def read_selfsm(path: str | Path) -> dict[str, float]:
    """``{sample: FREEMIX}`` from a VerifyBamID .selfSM file."""
    return {row["#SEQ_ID"]: float(row["FREEMIX"]) for row in _table(path)}


def read_sexcheck(path: str | Path) -> dict[str, tuple[int, float]]:
    """``{IID: (reported sex, X-chromosome F)}`` from a PLINK .sexcheck file.

    Rows where PLINK could not estimate F (``nan``) are skipped.
    """
    out: dict[str, tuple[int, float]] = {}
    for row in _table(path):
        f = float(row["F"])
        if f == f:  # skip NaN
            out[row["IID"]] = (int(row["PEDSEX"]), f)
    return out


def read_mosdepth_summaries(paths: Iterable[str | Path]) -> dict[str, float]:
    """``{sample: mean depth}`` from mosdepth summary files.

    The sample name is the file name before ``.mosdepth.summary.txt``.
    """
    out: dict[str, float] = {}
    for path in paths:
        name = Path(path).name.removesuffix(".mosdepth.summary.txt")
        for row in _table(path):
            if row["chrom"] == "total":
                out[name] = float(row["mean"])
    return out


def read_lmiss(path: str | Path) -> dict[str, float]:
    """``{SNP: call rate}`` from a PLINK .lmiss file."""
    return {row["SNP"]: 1.0 - float(row["F_MISS"]) for row in _table(path)}
