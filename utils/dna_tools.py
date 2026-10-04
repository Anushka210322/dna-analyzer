"""Core DNA analysis functions.

All public functions are case-insensitive and ignore surrounding whitespace.
Only the four standard DNA bases (A, C, G, T) are considered valid.
"""

VALID_BASES = frozenset("ACGT")

COMPLEMENT_MAP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}


def _clean(sequence):
    """Strip whitespace and convert a sequence to uppercase.

    Args:
        sequence (str): The raw sequence.

    Returns:
        str: The cleaned, uppercase sequence with internal whitespace removed.

    Raises:
        TypeError: If ``sequence`` is not a string.
    """
    if not isinstance(sequence, str):
        raise TypeError(
            f"Sequence must be a string, got {type(sequence).__name__}."
        )
    return "".join(sequence.split()).upper()


def validate_dna(sequence):
    """Check whether a sequence contains only valid DNA bases (A, C, G, T).

    An empty sequence (or whitespace-only input) is considered invalid.

    Args:
        sequence (str): The DNA sequence to check.

    Returns:
        bool: True if the sequence is valid, False otherwise.

    Examples:
        >>> validate_dna("ATGC")
        True
        >>> validate_dna("atgc")
        True
        >>> validate_dna("ATGX")
        False
    """
    if not isinstance(sequence, str):
        return False

    cleaned = _clean(sequence)
    if not cleaned:
        return False

    return set(cleaned) <= VALID_BASES


def _require_valid_dna(sequence):
    """Return a cleaned sequence or raise ValueError if it is invalid."""
    if not validate_dna(sequence):
        raise ValueError(f"Invalid DNA sequence: {sequence!r}")
    return _clean(sequence)


def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence as a percentage.

    GC content is the proportion of bases that are either G or C.

    Args:
        sequence (str): A valid DNA sequence.

    Returns:
        float: GC content as a percentage between 0.0 and 100.0.

    Raises:
        ValueError: If the sequence is not valid DNA.

    Examples:
        >>> calculate_gc_content("GGCC")
        100.0
        >>> calculate_gc_content("ATGC")
        50.0
    """
    cleaned = _require_valid_dna(sequence)
    gc_count = cleaned.count("G") + cleaned.count("C")
    return (gc_count / len(cleaned)) * 100


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence.

    The complement swaps A<->T and C<->G, and the result is then reversed,
    which gives the sequence of the opposite DNA strand read 5' to 3'.

    Args:
        sequence (str): A valid DNA sequence.

    Returns:
        str: The reverse complement, in uppercase.

    Raises:
        ValueError: If the sequence is not valid DNA.

    Examples:
        >>> reverse_complement("ATGC")
        'GCAT'
        >>> reverse_complement("AAAA")
        'TTTT'
    """
    cleaned = _require_valid_dna(sequence)
    complement = "".join(COMPLEMENT_MAP[base] for base in cleaned)
    return complement[::-1]


def calculate_base_counts(sequence):
    """Count occurrences of each DNA base in a sequence.

    Args:
        sequence (str): A valid DNA sequence.

    Returns:
        dict: Mapping of base letters ``A``, ``C``, ``G``, ``T`` to integer counts.
            Missing bases are included with a count of ``0``.

    Raises:
        ValueError: If the sequence is not valid DNA.

    Examples:
        >>> calculate_base_counts("ATGC")
        {'A': 1, 'C': 1, 'G': 1, 'T': 1}
        >>> calculate_base_counts("AAA")
        {'A': 3, 'C': 0, 'G': 0, 'T': 0}
    """
    cleaned = _require_valid_dna(sequence)
    return {base: cleaned.count(base) for base in "ACGT"}
