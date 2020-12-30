#!/usr/bin/env python3
"""Interfaces for reading and writing data also known as input/output (I/O)
"""
import Bio.SeqIO
import gzip


VALID_FORMATS = [
    "fasta",
    "genbank"
]


def infer_sequence_file_type(path_or_buffer):
    """Infer the file type of sequences in the given path or buffer.

    Parameters
    ----------
    path_or_buffer : str or handle
       A path to or handle for a sequence file whose type is supported by BioPython


    Raises
    ------
    ValueError
        If the given path's format cannot be determined

    Returns
    -------
    str :
        format type supported by BioPython's SeqIO parse interface
    """
    path_format = None
    for format in VALID_FORMATS:
        try:
            # Try to read the first record in the given path to determine its
            # file format. Work from the most to least common formats.
            sequence = next(Bio.SeqIO.parse(path_or_buffer, format))
            path_format = format
            break
        except StopIteration:
            pass

    if path_format is None:
        raise ValueError(f"Could not determine format for sequence path: '{path}'")

    try:
        path_or_buffer.seek(0)
    except AttributeError:
        pass

    return path_format


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
        compression = False

        try:
            path_format = infer_sequence_file_type(path)
        except UnicodeDecodeError:
            compression = True
            handle = gzip.open(path, "rt")
            path_format = infer_sequence_file_type(handle)
            path = handle

        sequences = Bio.SeqIO.parse(path, path_format)

        for sequence in sequences:
            yield sequence

        if compression:
            handle.close()
