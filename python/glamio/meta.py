"""Reader for the plain-text ``key = value`` sidecar files written by GLAM.

A sidecar (``*.meta``) holds one ``key = value`` pair per line.  The line is
split on the first ``=``; both sides are stripped; values stay strings.  Blank
lines and lines without ``=`` are ignored.  Values may carry a trailing comment
after the number, e.g. ``a_v = 0.4973  (velocity epoch)``, which is why
:func:`meta_float` converts only the first whitespace token.
"""

from __future__ import annotations

import os


def read_meta(path: str | os.PathLike) -> dict[str, str]:
    """Parse a ``key = value`` sidecar into a ``{key: value}`` dict of strings.

    Lines are split on the FIRST ``=``; keys and values are stripped.  Blank
    lines, lines without ``=`` and lines with an empty key are skipped.  A
    repeated key keeps the last value.
    """
    out: dict[str, str] = {}
    with open(os.fspath(path), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            if not key:
                continue
            out[key] = value.strip()
    return out


def meta_float(meta: dict[str, str], key: str, default: float | None = None) -> float | None:
    """Return ``float(first token of meta[key])``, or ``default`` if the key is absent.

    Only the first whitespace-separated token is converted, so a value such as
    ``"0.4973  velocity epoch"`` gives ``0.4973``.  Fortran ``D`` exponents
    (``1.0D-5``) are accepted.  An absent key or an empty value yields
    ``default``; a present but non-numeric first token raises ``ValueError``.
    """
    value = meta.get(key)
    if value is None:
        return default
    tokens = value.split()
    if not tokens:
        return default
    token = tokens[0].replace("D", "E").replace("d", "e")
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(f"meta key {key!r}: cannot convert {tokens[0]!r} to float") from exc
