"""Helpers for multiple sequence alignments held as lists of strings."""

from collections import Counter
from typing import Sequence


def all_strings_same_length(strings: Sequence[str]) -> bool:
    """True if every string has the same length (an empty list returns False)."""
    return len(set(map(len, strings))) == 1


def consensus_sequence(msa: Sequence[str], threshold: float, keep_gaps: bool = False,
                       min_count: int = 0) -> str:
    """Majority-rule consensus of an alignment, ignoring gaps (``-``) when counting.

    For each column the most common non-gap character is used if
    (its frequency among non-gap characters) >= ``threshold`` and its count is
    >= ``min_count``.

    Parameters
    ----------
    msa : list of str
        Aligned sequences of equal length.
    threshold : float
        Minimum fraction (0-1) of the most common character.
    keep_gaps : bool
        What to do with a column that fails the test (or contains only gaps):
        ``True`` writes ``-`` so the consensus keeps the alignment length;
        ``False`` (default) drops the column, so the consensus can be shorter.
    min_count : int
        Minimum number of sequences supporting the consensus character.
    """
    consensus = []
    for column in zip(*msa):
        counter = Counter(column)
        del counter["-"]
        seq_number = sum(counter.values())
        if seq_number:
            char, count = counter.most_common(1)[0]
            if count / seq_number >= threshold and count >= min_count:
                consensus.append(char)
                continue
        if keep_gaps:
            consensus.append("-")
    return "".join(consensus)


def consensus_sequence_with_dashes(msa: Sequence[str], threshold: float, min_count: int = 5) -> str:
    """Consensus that keeps the alignment length (``-`` for unsupported columns).

    Equivalent to ``consensus_sequence(msa, threshold, keep_gaps=True, min_count=5)``.
    """
    return consensus_sequence(msa, threshold, keep_gaps=True, min_count=min_count)
