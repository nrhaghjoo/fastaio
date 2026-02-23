import os
from typing import Literal

def read_fasta(file_path, seq_type = Literal["DNA", "RNA", "protein"], is_alignment = False):
    """Reads a FASTA file and validates all possible errors before processing."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")

    headers = []
    sequences = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"Error: The file '{file_path}' is empty.")
        show_error(f"Error: The file '{file_path}' is empty.")

    if seq_type == "protein":
        valid_chars = "ACDEFGHIKLMNPQRSTVWY"
    else:
        valid_chars = "TCGAU"

    def all_strings_same_length(strings):
        return len(set(map(len, strings))) == 1

    if is_alignment:
        if not all_strings_same_length(sequences):
            raise ValueError(f"Error: The file '{file_path}' is not a multiple sequence alignment.")


    sequence_data = ''
    for line in lines:
        line = line.strip()

        if line.startswith(">"):
            if sequence_data:
                sequence_data = re.sub(f"[^{valid_chars}]", "-", sequence_data.upper()).replace("U", "T")
                sequences.append(sequence_data)
                sequence_data = ''

            headers.append(line[1:])  # Store header (without '>')

        elif line:
            sequence_data += line


    if sequence_data:
        sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
        sequences.append(sequence_data)

    if not headers:
        raise ValueError(f"Error: The file '{file_path}' does not contain any headers.")

    if len(headers) != len(sequences):
        raise ValueError(f"Error: Mismatch between headers ({len(headers)}) and sequences ({len(sequences)}).")


    return headers, sequences