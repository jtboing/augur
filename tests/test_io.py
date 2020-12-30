#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import lzma
from pathlib import Path
import pytest
import random
import sys

# make sure we get the local version of modules
sys.path.insert(0, (str(Path(__file__).parent.parent.parent)))
from augur.io import read_sequences


@pytest.fixture
def sequences():
    def random_seq(k):
        return "".join(random.choices(("A","T","G","C"), k=k))
    return [
        SeqRecord(Seq(random_seq(10)), id="SEQ_1"),
        SeqRecord(Seq(random_seq(10)), id="SEQ_2"),
        SeqRecord(Seq(random_seq(10)), id="SEQ_3")
    ]

@pytest.fixture
def fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta")
    SeqIO.write(sequences, filename, "fasta")
    return filename

@pytest.fixture
def additional_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "additional_sequences.fasta")
    SeqIO.write(sequences, filename, "fasta")
    return filename

@pytest.fixture
def gzip_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta.gz")

    with gzip.open(filename, "wt") as oh:
        SeqIO.write(sequences, oh, "fasta")

    return filename

@pytest.fixture
def lzma_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta.xz")

    with lzma.open(filename, "wt") as oh:
        SeqIO.write(sequences, oh, "fasta")

    return filename

@pytest.fixture
def genbank_reference():
    return "tests/builds/zika/config/zika_outgroup.gb"


class TestReadSequences:
    def test_read_sequences_from_single_file(self, fasta_filename):
        sequences = read_sequences(fasta_filename, "fasta")
        assert len([sequence for sequence in sequences]) == 3

    def test_read_sequences_from_multiple_files(self, fasta_filename, additional_fasta_filename):
        sequences = read_sequences([fasta_filename, additional_fasta_filename], "fasta")
        assert len([sequence for sequence in sequences]) == 6

    def test_read_single_fasta_record(self, fasta_filename):
        record = next(read_sequences(fasta_filename, "fasta"))
        assert record.id == "SEQ_1"

    def test_read_single_genbank_record(self, genbank_reference):
        reference = next(read_sequences(genbank_reference, "genbank"))
        assert reference.id == "KX369547.1"

    def test_read_single_genbank_record_from_a_path(self, genbank_reference):
        reference = next(read_sequences(Path(genbank_reference), "genbank"))
        assert reference.id == "KX369547.1"

    def test_read_sequences_from_single_gzip_file(self, gzip_fasta_filename):
        sequences = read_sequences(gzip_fasta_filename, "fasta")
        assert len([sequence for sequence in sequences]) == 3

    def test_read_sequences_from_single_lzma_file(self, lzma_fasta_filename):
        sequences = read_sequences(lzma_fasta_filename, "fasta")
        assert len([sequence for sequence in sequences]) == 3

    def test_read_sequences_from_multiple_files_with_different_compression(self, fasta_filename, gzip_fasta_filename, lzma_fasta_filename):
        sequences = read_sequences([fasta_filename, gzip_fasta_filename, lzma_fasta_filename], "fasta")
        assert len([sequence for sequence in sequences]) == 9
