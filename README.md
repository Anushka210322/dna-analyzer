# DNA Sequence Analyzer

A modular Python command-line tool for basic DNA sequence analysis.
It validates a sequence, counts bases, calculates GC content, and computes
the reverse complement. Input can be typed on the command line or read from
a FASTA or plain-text file.

## Features

- **Validation**: confirms a sequence contains only A, C, G, T (case-insensitive)
- **GC content**: percentage of G and C bases
- **Reverse complement**: sequence of the opposite DNA strand (5' to 3')
- **Base counts**: number of A, C, G, and T residues
- **File support**: plain-text files and multi-record FASTA files
- **Verbose output** (`-v`): length, validity, and per-base percentages
- **File output** (`-o`): write the report to a path you choose
- **Unit tests** written with `pytest`

## Project Structure

```
dna-analyzer/
├── .gitignore
├── requirements.txt
├── README.md
├── main.py                  # CLI entry point (argparse)
├── utils/
│   ├── __init__.py
│   └── dna_tools.py         # Core analysis functions
└── tests/
    ├── __init__.py
    └── test_dna_tools.py    # pytest unit tests
```

## Installation

Requires **Python 3.8 or newer**.

1. Open a terminal in the project folder:

```bash
cd dna-analyzer
```

2. (Recommended) Create and activate a virtual environment:

```bash
python -m venv venv

# macOS / Linux
source venv/bin/activate

# Windows (PowerShell)
venv\Scripts\Activate.ps1
```

3. Install the dependencies:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python main.py (-s SEQUENCE | -f FILE) [--gc] [--revcomp] [--all] [-v] [-o FILE]
```

| Flag | Description |
|------|-------------|
| `-s`, `--sequence` | DNA sequence typed directly on the command line |
| `-f`, `--file` | Path to a FASTA or plain-text file |
| `--gc` | Calculate GC content |
| `--revcomp` | Compute the reverse complement |
| `--all` | Run all analyses (validation, GC, reverse complement, base counts) |
| `-v`, `--verbose` | Extra details (length, validity, per-base percentages) |
| `-o`, `--output` | Write the report to a file as well as printing it |

You must provide **either** `-s` **or** `-f` (not both), plus at least one of
`--gc`, `--revcomp`, or `--all`.

### Examples

GC content of a sequence:

```bash
python main.py -s ATGCGC --gc
```

```
Sequence: command-line input
  GC content: 66.67%
```

Reverse complement:

```bash
python main.py -s ATGCGC --revcomp
```

Everything at once, with verbose details, writing a report file:

```bash
python main.py -s ATGCGC --all -v -o report.txt
```

From a FASTA file:

```bash
python main.py -f sequences.fasta --all
```

Example `sequences.fasta`:

```
>gene_1
ATGCGTACGTTAGC
>gene_2
GGCCGGCCAATT
```

Show help:

```bash
python main.py -h
```

### Exit codes

- `0`: all sequences were valid and analyzed
- `1`: a sequence was invalid, no analysis flag was given, or a file could not be read/written

## Running the Tests

From the project root:

```bash
pytest
```

For more detail:

```bash
pytest -v
```

## How It Works

- **GC content** = (number of G + number of C) / sequence length × 100
- **Reverse complement**: each base is swapped for its partner
  (A↔T, C↔G) and the result is reversed
- **Base counts**: occurrences of A, C, G, and T after the sequence is
  uppercased and whitespace is removed

## License

Free to use and modify for learning and personal projects.
