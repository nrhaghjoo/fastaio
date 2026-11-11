import os

import pytest

from fastaio import (PATH, all_strings_same_length, consensus_sequence, consensus_sequence_with_dashes,
                     df_2_fasta, read_fasta, read_reference, save_2_fasta, separate_fasta)


def write(tmp_path, text, name="x.fasta"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_dna_cleaning_and_multiline(tmp_path):
    f = write(tmp_path, ">s1 desc\nacgu\nNN\n>s2\nACGTRY\n")
    h, s = read_fasta(f)
    assert h == ["s1 desc", "s2"]
    assert s == ["ACGT--", "ACGT--"]          # U->T, ambiguity -> '-'


def test_last_record_uses_same_rule_as_others(tmp_path):
    f = write(tmp_path, ">p1\nmkv-x\n>p2\nmkl-y\n")
    _, s = read_fasta(f, seq_type="protein")
    assert s == ["MKV-X", "MKL-Y"]


def test_raw_keeps_case(tmp_path):
    f = write(tmp_path, ">a\nAcG\n>b\nAtG\n")
    assert read_fasta(f, seq_type="raw")[1] == ["AcG", "AtG"]


def test_alignment_check(tmp_path):
    ok = write(tmp_path, ">a\nACGT\n>b\nAC-T\n", "ok.fasta")
    bad = write(tmp_path, ">a\nACGT\n>b\nACG\n", "bad.fasta")
    assert read_fasta(ok, is_alignment=True)[1] == ["ACGT", "AC-T"]
    read_fasta(bad)                                   # plain FASTA: fine
    with pytest.raises(ValueError, match="not a multiple sequence alignment"):
        read_fasta(bad, is_alignment=True)


def test_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_fasta(tmp_path / "missing.fasta")
    with pytest.raises(ValueError, match="empty"):
        read_fasta(write(tmp_path, "", "e.fasta"))
    with pytest.raises(ValueError, match="headers"):
        read_fasta(write(tmp_path, "ACGT\n", "n.fasta"))
    with pytest.raises(ValueError, match="seq_type"):
        read_fasta(write(tmp_path, ">a\nA\n", "t.fasta"), seq_type="dna")


def test_read_reference(tmp_path):
    assert read_reference(write(tmp_path, ">ref\nACGT\n>q\nTTTT\n")) == "ACGT"


def test_consensus_variants():
    msa = ["AC-T", "AC-T", "AG-T", "TC-A"]
    assert consensus_sequence(msa, 0.6) == "ACT"                     # drops gap column
    assert consensus_sequence(msa, 0.6, keep_gaps=True) == "AC-T"
    assert consensus_sequence(msa, 0.6, keep_gaps=True, min_count=4) == "----"
    assert consensus_sequence_with_dashes(["A"] * 5, 0.5) == "A"
    assert consensus_sequence_with_dashes(["A"] * 4, 0.5) == "-"     # fewer than 5 supporting
    assert all_strings_same_length(msa) and not all_strings_same_length([])


def test_writers_roundtrip(tmp_path):
    pd = pytest.importorskip("pandas")
    out = save_2_fasta(["a", "b"], ["AC", "GT"], str(tmp_path / "o.fasta"), verbose=False)
    assert read_fasta(out) == (["a", "b"], ["AC", "GT"])
    df = pd.DataFrame({"id": ["a"], "seq": ["AC"]})
    assert read_fasta(df_2_fasta(df, str(tmp_path / "d.fasta"), "id", "seq"))[1] == ["AC"]
    paths = separate_fasta([0, 1, 0, -1], list("wxyz"), ["A", "C", "G", "T"], str(tmp_path / "cl"), "s.fasta")
    assert sorted(os.path.basename(p) for p in paths) == ["cluster_-1_s.fasta", "cluster_0_s.fasta", "cluster_1_s.fasta"]
    assert read_fasta(str(tmp_path / "cl" / "cluster_0_s.fasta"))[0] == ["w", "y"]
    with pytest.raises(ValueError):
        save_2_fasta(["a"], [], str(tmp_path / "z.fasta"))


def test_path_helper(tmp_path):
    p = PATH(str(tmp_path / "data" / "sample1.fasta"))
    assert p.add_to_name("distance", extention=".csv") == str(tmp_path / "data" / "sample1_distance.csv")
    p.go_up(1).add_to_dir("result").create_folder()
    assert os.path.isdir(tmp_path / "result")
    assert PATH("s.fasta").add_to_dir("result").dir == "result"     # no stray leading '/'
