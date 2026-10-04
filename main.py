#!/usr/bin/env python3
"""Command-line interface for DNA sequence analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from utils.dna_tools import (
    calculate_base_counts,
    calculate_gc_content,
    reverse_complement,
    validate_dna,
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Analyze DNA sequences: validate bases, GC content, "
        "reverse complement, and base counts.",
        epilog="Provide either -s/--sequence or -f/--file, plus at least one "
        "of --gc, --revcomp, or --all.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-s",
        "--sequence",
        help="DNA sequence typed on the command line",
    )
    source.add_argument(
        "-f",
        "--file",
        type=Path,
        help="Path to a FASTA or plain-text file of DNA sequences",
    )
    parser.add_argument(
        "--gc",
        action="store_true",
        help="Calculate GC content",
    )
    parser.add_argument(
        "--revcomp",
        action="store_true",
        help="Compute the reverse complement",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all analyses (validation, GC content, reverse complement, "
        "base counts)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print extra details (length, per-base percentages, file path)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write the report to this file in addition to printing it",
    )
    return parser.parse_args(argv)


def read_sequences_from_file(path):
    """Load named sequences from a FASTA or plain-text file.

    FASTA records start with a ``>`` header line. Plain-text files are treated
    as a single sequence (whitespace is ignored later by the analysis functions).

    Returns:
        list[tuple[str, str]]: (label, raw_sequence) pairs.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is empty or contains no sequence data.
    """
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    text = path.read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"File is empty: {path}")

    lines = text.splitlines()
    fasta_headers = [line for line in lines if line.startswith(">")]
    if fasta_headers:
        records = []
        current_name = None
        current_chunks = []
        for line in lines:
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, "".join(current_chunks)))
                current_name = line[1:].strip() or f"record_{len(records) + 1}"
                current_chunks = []
            else:
                current_chunks.append(line.strip())
        if current_name is not None:
            records.append((current_name, "".join(current_chunks)))
        nonempty = [(name, seq) for name, seq in records if "".join(seq.split())]
        if not nonempty:
            raise ValueError(f"No sequence data found in FASTA file: {path}")
        return nonempty

    return [(path.name, text)]


def _normalized(sequence):
    """Uppercase DNA and drop all whitespace."""
    if not isinstance(sequence, str):
        return ""
    return "".join(sequence.split()).upper()


def format_report(label, sequence, show_gc, show_revcomp, show_all, verbose):
    """Build a human-readable analysis report for one sequence."""
    lines = [f"Sequence: {label}"]
    is_valid = validate_dna(sequence)
    cleaned = _normalized(sequence)

    if verbose or show_all:
        lines.append(f"  Valid DNA : {'Yes' if is_valid else 'No'}")
        lines.append(f"  Length    : {len(cleaned)} bases")

    if not is_valid:
        lines.append("  Error     : invalid DNA (only A, C, G, T are allowed)")
        return "\n".join(lines), False

    if show_gc or show_all:
        gc = calculate_gc_content(sequence)
        lines.append(f"  GC content: {gc:.2f}%")

    if show_revcomp or show_all:
        rc = reverse_complement(sequence)
        lines.append(f"  Reverse complement: {rc}")

    if show_all:
        counts = calculate_base_counts(sequence)
        count_text = ", ".join(f"{base}={counts[base]}" for base in "ACGT")
        lines.append(f"  Base counts: {count_text}")
        if verbose:
            total = len(cleaned)
            pct = ", ".join(
                f"{base}={counts[base] / total * 100:.2f}%" for base in "ACGT"
            )
            lines.append(f"  Base %    : {pct}")

    return "\n".join(lines), True


def run(argv=None):
    args = parse_args(argv)

    if not (args.gc or args.revcomp or args.all):
        print(
            "error: specify at least one of --gc, --revcomp, or --all",
            file=sys.stderr,
        )
        return 1

    try:
        if args.sequence is not None:
            records = [("command-line input", args.sequence)]
        else:
            records = read_sequences_from_file(args.file)
            if args.verbose:
                print(f"Reading sequences from: {args.file.resolve()}")
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    reports = []
    all_ok = True
    for label, sequence in records:
        report, ok = format_report(
            label,
            sequence,
            show_gc=args.gc,
            show_revcomp=args.revcomp,
            show_all=args.all,
            verbose=args.verbose,
        )
        reports.append(report)
        all_ok = all_ok and ok

    output_text = "\n\n".join(reports) + "\n"
    print(output_text, end="")

    if args.output is not None:
        try:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(output_text, encoding="utf-8")
            if args.verbose:
                print(f"Wrote report to: {args.output.resolve()}")
        except OSError as exc:
            print(f"error: could not write output file: {exc}", file=sys.stderr)
            return 1

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(run())
