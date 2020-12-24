#!/usr/bin/env python3
"""Interfaces for reading and writing data also known as input/output (I/O)
"""
import Bio.SeqIO


def read_sequences(paths_or_buffers):
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc.) and file types
    (e.g., FASTA, GenBank, etc.) and return a stream of sequence records.

    Parameters
    ----------
    paths_or_buffers : str, handle, or list of str or handle
        One or more paths to sequence files of any type supported by BioPython.

    Yields
    ------
    Bio.SeqRecord.SeqRecord
        Sequence record from the given path(s).
    """
    if isinstance(paths_or_buffers, str) or hasattr(paths_or_buffers, "read"):
        paths_or_buffers = [paths_or_buffers]

    for path in paths_or_buffers:
        # TODO: infer file type for the given path.
        sequences = Bio.SeqIO.parse(path, "fasta")

        for sequence in sequences:
            yield sequence
