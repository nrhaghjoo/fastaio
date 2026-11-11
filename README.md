# fastaio

Small, dependency-free helpers for FASTA files and multiple sequence alignments
(MSAs), shared by my HIV/SIV mutation-analysis projects.

## Install

```bash
pip install git+https://github.com/nrhaghjoo/fastaio.git
```

> The name `fastaio` on PyPI belongs to an unrelated project, so install from GitHub as above.

## Reading FASTA files

```python
from fastaio import read_fasta, read_reference

# nucleotide alignment: every sequence must have the same length
headers, seqs = read_fasta("reads_aligned.fasta", seq_type="DNA", is_alignment=True)

# plain (unaligned) nucleotide FASTA
headers, seqs = read_fasta("contigs.fasta", seq_type="DNA")

# protein alignment, upper-cased, gaps kept
headers, seqs = read_fasta("family.fasta", seq_type="protein", is_alignment=True)

# first sequence only (e.g. a reference genome)
ref = read_reference("HXB2.fasta")
```

| `seq_type` | What happens to each sequence |
|---|---|
| `"DNA"` / `"RNA"` (default) | upper-case; any character other than A, C, G, T, U becomes `-`; U → T |
| `"protein"` | upper-case only |
| `"raw"` | unchanged (keeps lower-case insert states) |

`is_alignment=True` raises `ValueError` when sequences differ in length.

## Writing

```python
from fastaio import save_2_fasta, df_2_fasta, separate_fasta

save_2_fasta(headers, seqs, "out.fasta")
df_2_fasta(df, "out.fasta", header="contig", sequence="sequences")   # needs pandas
separate_fasta(labels, headers, seqs, "clusters/", "sample.fasta")   # one file per cluster label
```

## Alignment helpers

```python
from fastaio import consensus_sequence, consensus_sequence_with_dashes

consensus_sequence(seqs, threshold=0.5)                       # drops unsupported columns
consensus_sequence(seqs, 0.5, keep_gaps=True, min_count=5)    # keeps alignment length
consensus_sequence_with_dashes(seqs, 0.5)                     # same as the line above
```

## Output paths

```python
from fastaio import PATH

p = PATH("data/sample1.fasta")
p.add_to_name("distance", extention=".csv")        # data/sample1_distance.csv
p.go_up(1).add_to_dir("result").create_folder()    # creates ./result
```

## Tests

```bash
pip install -e ".[test]"
pytest
```
