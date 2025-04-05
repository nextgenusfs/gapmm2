#!/usr/bin/env python3

import pytest
from gapmm2.align import cs2tuples, cs2coords


def test_cs2tuples():
    """Test cs2tuples function with various CIGAR strings."""
    # Test with a simple CIGAR string
    cs = ":100"
    result = list(cs2tuples(cs))
    assert result == [(":", "100")]

    # Test with a more complex CIGAR string
    cs = ":50*at:20+aaa-ggg:30~gt100ag:40"
    result = list(cs2tuples(cs))
    assert result == [
        (":", "50"),
        ("*", "at"),
        (":", "20"),
        ("+", "aaa"),
        ("-", "ggg"),
        (":", "30"),
        ("~", "gt100ag"),
        (":", "40"),
    ]

    # Test with an empty string
    cs = ""
    result = list(cs2tuples(cs))
    assert result == []


def test_cs2coords_forward_strand():
    """Test cs2coords function with forward strand alignments."""
    # Simple case with one exon
    start = 1000
    qstart = 0
    length = 100
    strand = "+"
    cs = ":100"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1001, 1100)]  # 1-based coordinates with offset=1
    assert queries == [(0, 100)]
    assert mismatches == 0
    assert gaps == 0
    assert proper_splice is True

    # Case with two exons and a proper splice junction
    start = 1000
    qstart = 0
    length = 200
    strand = "+"
    cs = ":100~gt50ag:100"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1001, 1100), (1151, 1250)]
    assert queries == [(0, 100), (101, 200)]
    assert mismatches == 0
    assert gaps == 0
    assert proper_splice is True

    # Case with a mismatch
    start = 1000
    qstart = 0
    length = 100
    strand = "+"
    cs = ":50*at:49"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1001, 1099)]
    assert queries == [(0, 100)]
    assert mismatches == 1
    assert gaps == 0
    assert proper_splice is True


def test_cs2coords_reverse_strand():
    """Test cs2coords function with reverse strand alignments."""
    # Simple case with one exon
    start = 1000
    qstart = 0
    length = 100
    strand = "-"
    cs = ":100"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1001, 1100)]
    assert queries == [(0, 100)]
    assert mismatches == 0
    assert gaps == 0
    assert proper_splice is True

    # Case with two exons and a proper splice junction
    start = 1000
    qstart = 0
    length = 200
    strand = "-"
    cs = ":100~ct50ac:100"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1151, 1250), (1001, 1100)]
    assert queries == [(0, 100), (101, 200)]
    assert mismatches == 0
    assert gaps == 0
    assert proper_splice is True


def test_cs2coords_with_indels():
    """Test cs2coords function with insertions and deletions."""
    # Case with an insertion
    start = 1000
    qstart = 0
    length = 103
    strand = "+"
    cs = ":50+aaa:50"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1001, 1103)]
    assert queries == [(0, 103)]
    assert mismatches == 0
    assert gaps == 1
    assert proper_splice is True

    # Case with a deletion
    start = 1000
    qstart = 0
    length = 97
    strand = "+"
    cs = ":50-ggg:50"

    exons, queries, mismatches, gaps, proper_splice = cs2coords(
        start, qstart, length, strand, cs
    )

    assert exons == [(1001, 1103)]
    assert queries == [(0, 97)]
    assert mismatches == 0
    assert gaps == 1
    assert proper_splice is True
