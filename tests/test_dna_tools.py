"""Unit tests for utils/dna_tools.py (run with: pytest)."""

import pytest

from utils.dna_tools import (
    calculate_base_counts,
    calculate_gc_content,
    reverse_complement,
    validate_dna,
)


class TestValidateDna:
    def test_valid_uppercase(self):
        assert validate_dna("ATGC") is True

    def test_valid_lowercase(self):
        assert validate_dna("atgc") is True

    def test_valid_mixed_case(self):
        assert validate_dna("AtGc") is True

    def test_valid_with_whitespace(self):
        assert validate_dna("  ATG C\n") is True

    def test_invalid_character(self):
        assert validate_dna("ATGX") is False

    def test_rna_base_is_invalid(self):
        assert validate_dna("AUGC") is False

    def test_digits_are_invalid(self):
        assert validate_dna("ATG1") is False

    def test_empty_string_is_invalid(self):
        assert validate_dna("") is False

    def test_whitespace_only_is_invalid(self):
        assert validate_dna("   ") is False

    def test_non_string_is_invalid(self):
        assert validate_dna(None) is False
        assert validate_dna(123) is False
        assert validate_dna(["A", "T"]) is False


class TestCalculateGcContent:
    def test_all_gc(self):
        assert calculate_gc_content("GGCC") == 100.0

    def test_no_gc(self):
        assert calculate_gc_content("ATAT") == 0.0

    def test_half_gc(self):
        assert calculate_gc_content("ATGC") == 50.0

    def test_lowercase_input(self):
        assert calculate_gc_content("atgc") == 50.0

    def test_whitespace_is_ignored(self):
        assert calculate_gc_content("AT GC") == 50.0

    def test_fractional_result(self):
        assert calculate_gc_content("GAT") == pytest.approx(33.3333, abs=1e-3)

    def test_single_base_g(self):
        assert calculate_gc_content("G") == 100.0

    def test_single_base_a(self):
        assert calculate_gc_content("A") == 0.0

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            calculate_gc_content("ATGX")

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError):
            calculate_gc_content("")

    def test_non_string_raises(self):
        with pytest.raises(ValueError):
            calculate_gc_content(None)


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_all_a(self):
        assert reverse_complement("AAAA") == "TTTT"

    def test_all_complements(self):
        assert reverse_complement("ACGT") == "ACGT"

    def test_palindromic_sequence(self):
        # GAATTC (EcoRI site) is its own reverse complement.
        assert reverse_complement("GAATTC") == "GAATTC"

    def test_lowercase_input_returns_uppercase(self):
        assert reverse_complement("atgc") == "GCAT"

    def test_whitespace_is_ignored(self):
        assert reverse_complement("AT GC") == "GCAT"

    def test_double_reverse_complement_returns_original(self):
        original = "ATGCGTACCGTA"
        assert reverse_complement(reverse_complement(original)) == original

    def test_single_base(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            reverse_complement("ATGX")

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError):
            reverse_complement("")


class TestCalculateBaseCounts:
    def test_balanced(self):
        assert calculate_base_counts("ATGC") == {"A": 1, "C": 1, "G": 1, "T": 1}

    def test_all_one_base(self):
        assert calculate_base_counts("AAA") == {"A": 3, "C": 0, "G": 0, "T": 0}

    def test_missing_bases_are_zero(self):
        counts = calculate_base_counts("GGCC")
        assert counts["A"] == 0
        assert counts["T"] == 0
        assert counts["G"] == 2
        assert counts["C"] == 2

    def test_lowercase_input(self):
        assert calculate_base_counts("atgc") == {"A": 1, "C": 1, "G": 1, "T": 1}

    def test_whitespace_is_ignored(self):
        assert calculate_base_counts("AA TT") == {"A": 2, "C": 0, "G": 0, "T": 2}

    def test_single_base(self):
        assert calculate_base_counts("C") == {"A": 0, "C": 1, "G": 0, "T": 0}

    def test_keys_are_always_acgt(self):
        counts = calculate_base_counts("T")
        assert list(counts.keys()) == ["A", "C", "G", "T"]

    def test_invalid_letter_raises(self):
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            calculate_base_counts("ATGX")

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError):
            calculate_base_counts("")

    def test_rna_u_raises(self):
        with pytest.raises(ValueError):
            calculate_base_counts("AUGC")
