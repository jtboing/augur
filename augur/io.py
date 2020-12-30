#!/usr/bin/env python3
"""Interfaces for reading and writing data also known as input/output (I/O)
"""
import Bio.SeqIO
import Bio.SeqRecord
from pathlib import Path
from xopen import xopen


def read_sequences(paths, format):
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc.) and return a stream
    of sequence records in the requested format (e.g., "fasta", "genbank", etc.).

    Parameters
    ----------
    paths : str, Path-like, or list of str or Path-like objects
        One or more paths to sequence files of any type supported by BioPython.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.).

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


def write_sequences(sequences, path, format, handle=None, return_handle=False):
    """Write sequences to a given path in the given format.

    Automatically infer compression mode (e.g., gzip, etc.) based on the path's
    filename extension.

    Parameters
    ----------
    sequences : iterable of Bio.SeqRecord.SeqRecord objects
        A list-like collection of sequences to write

    path : str or Path-like object
        A path to a file to write the given sequences in the given format.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.)

    handle : IO buffer, optional
        An existing file handle for the given path to use when writing sequence
        records one by one as in a for loop.

    return_handle : bool, optional
        Return the handle created for writing to the given path instead of
        closing it.

    Returns
    -------
    int :
        Number of sequences written out to the given path.

    IO buffer, optional :
        A file handle that was opened to write to the given path and was
        requested by the calling code with the `return_handle` argument.

    """
    # Don't allow writing individual sequence records unless the caller has
    # explicitly requested that the file handle remain open for subsequent
    # writes.
    if not return_handle and isinstance(sequences, Bio.SeqRecord.SeqRecord):
        raise TypeError("Input sequences must be an iterable of BioPython SeqRecord objects unless return_handle is True.")

    if handle is None:
        handle = xopen(path, "wt")

    # Bio.SeqIO supports writing to the same handle multiple times for specific
    # file formats. For the formats we use, this function call should work for
    # both a newly opened file handle or one that is provided by the caller.
    # For more details see:
    # https://github.com/biopython/biopython/blob/25f5152f4aeefe184a323db25694fbfe0593f0e2/Bio/SeqIO/__init__.py#L233-L251
    sequences_written = Bio.SeqIO.write(
        sequences,
        handle,
        format
    )

    if return_handle:
        return sequences_written, handle
    else:
        handle.close()
        return sequences_written
