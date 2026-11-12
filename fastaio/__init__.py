"""fastaio - small FASTA / alignment utilities shared by Niloo Haghjoo's analysis projects."""

from .alignment import all_strings_same_length, consensus_sequence, consensus_sequence_with_dashes
from .io import df_2_fasta, read_fasta, read_reference, save_2_fasta, separate_fasta
from .paths import PATH

__version__ = "0.2.0"

__all__ = [
    "read_fasta", "read_reference", "save_2_fasta", "df_2_fasta", "separate_fasta",
    "all_strings_same_length", "consensus_sequence", "consensus_sequence_with_dashes",
    "PATH",
]
