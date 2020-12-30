#!/usr/bin/env python3
"""Interfaces for reading and writing data also known as input/output (I/O)
"""
import Bio.SeqIO
from pathlib import Path
from xopen import xopen


def read_sequences(paths, format):
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc. and return a stream
    of sequence records in the requested format (e.g., "fasta", "genbank", etc.).

    Parameters
    ----------
    paths : str, Path-like, or list of str or Path-like objects
        One or more paths to sequence files of any type supported by BioPython.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.)

    Yields
    ------
    Bio.SeqRecord.SeqRecord
        Sequence record from the given path(s).

    """
    try:
        # Since we accept any Path-like string as an input, we try to cast the
        # given path to a Path instance and then place it in a list. Note that
        # this approach could break if we eventually support arbitrary URLs or
        # paths to S3 buckets, etc.
        paths = [Path(paths)]
    except TypeError:
        # If the given input is not a Path-like string, we assume it is an
        # iterable sequence already.
        pass

    for path in paths:
        # Open the given path as a handle, inferring the file's compression.
        # This way we can pass a handle to BioPython's SeqIO interface
        # regardless of the compression mode.
        with xopen(path, "r") as handle:
            sequences = Bio.SeqIO.parse(handle, format)

            for sequence in sequences:
                yield sequence
