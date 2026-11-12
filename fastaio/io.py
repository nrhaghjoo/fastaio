"""Reading and writing FASTA files.

All readers return plain Python lists so they can be used with any downstream
library. The cleaning rules are explicit through ``seq_type`` so the same
function covers nucleotide alignments, plain nucleotide FASTA files and
protein alignments.
"""

import os
import re
from typing import Literal, Sequence

from .alignment import all_strings_same_length

SeqType = Literal["DNA", "RNA", "protein", "raw"]

_NUCLEOTIDES = "TCGAU"


def _clean(sequence: str, seq_type: SeqType) -> str:
    """Normalise one sequence according to its type.

    * ``DNA`` / ``RNA``: upper-case, every character that is not T, C, G, A or U
      becomes ``-`` and U is converted to T.
    * ``protein``: upper-case only; gaps and ambiguity codes are kept.
    * ``raw``: returned unchanged (keeps lower-case insert states in alignments).
    """
    if seq_type in ("DNA", "RNA"):
        return re.sub(f"[^{_NUCLEOTIDES}]", "-", sequence.upper()).replace("U", "T")
    if seq_type == "protein":
        return sequence.upper()
    if seq_type == "raw":
        return sequence
    raise ValueError(f"seq_type must be 'DNA', 'RNA', 'protein' or 'raw', got {seq_type!r}")


def read_fasta(file_path, seq_type: SeqType = "DNA", is_alignment: bool = False):
    """Read a FASTA file and return ``(headers, sequences)``.

    Parameters
    ----------
    file_path : str or path-like
        FASTA file to read.
    seq_type : {"DNA", "RNA", "protein", "raw"}
        How sequences are cleaned (see ``_clean``). Default ``"DNA"``.
    is_alignment : bool
        If True the file must be a multiple sequence alignment: every sequence
        must have the same length, otherwise ``ValueError`` is raised.

    Returns
    -------
    headers : list of str
        Header lines without the leading ``>``.
    sequences : list of str
        Cleaned sequences, in file order.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")

    with open(file_path, "r") as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"Error: The file '{file_path}' is empty.")

    headers, sequences = [], []
    sequence_data = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if headers:
                sequences.append(_clean(sequence_data, seq_type))
            headers.append(line[1:])
            sequence_data = ""
        elif line:
            sequence_data += line
    if headers:
        sequences.append(_clean(sequence_data, seq_type))

    if not headers:
        raise ValueError(f"Error: The file '{file_path}' does not contain any headers.")

    if is_alignment and not all_strings_same_length(sequences):
        raise ValueError(f"Error: The file '{file_path}' is not a multiple sequence alignment.")

    return headers, sequences


def read_reference(file_path, seq_type: SeqType = "DNA") -> str:
    """Return only the first sequence of a FASTA file (e.g. a reference genome)."""
    _, sequences = read_fasta(file_path, seq_type=seq_type)
    return sequences[0]


def save_2_fasta(headers: Sequence[str], sequences: Sequence[str], file_path, verbose: bool = True):
    """Write headers and sequences to a FASTA file and return the file path."""
    if len(headers) != len(sequences):
        raise ValueError("Error: the length of headers and sequences are not the same")
    with open(file_path, "w") as file:
        for header, sequence in zip(headers, sequences):
            file.write(f">{header}\n{sequence}\n")
    if verbose:
        print(f"File saved into {file_path}")
    return file_path


def df_2_fasta(df, file_name, header, sequence):
    """Write a pandas DataFrame to FASTA using the ``header`` and ``sequence`` columns."""
    with open(file_name, "w") as file:
        for _, row in df.iterrows():
            file.write(">" + str(row[header]) + "\n")
            file.write(str(row[sequence]) + "\n")
    return file_name


def separate_fasta(labels, headers, sequences, output_dir, file_name):
    """Write one FASTA file per label: ``<output_dir>/cluster_<label>_<file_name>``.

    Typical use: split reads into clusters after DBSCAN (label ``-1`` = noise).
    Returns the list of written paths.
    """
    if not (len(labels) == len(headers) == len(sequences)):
        raise ValueError("labels, headers and sequences must have the same length")
    os.makedirs(output_dir, exist_ok=True)
    written = []
    for label in sorted(set(labels), key=str):
        idx = [i for i, lab in enumerate(labels) if lab == label]
        path = os.path.join(output_dir, f"cluster_{label}_{file_name}")
        save_2_fasta([headers[i] for i in idx], [sequences[i] for i in idx], path, verbose=False)
        written.append(path)
    return written
