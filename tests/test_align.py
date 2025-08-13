#!/usr/bin/env python3

from unittest.mock import MagicMock, patch

from gapmm2.align import cs2coords, cs2tuples, splice_aligner, splice_aligner_minimap2


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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

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

    exons, queries, mismatches, gaps, proper_splice = cs2coords(start, qstart, length, strand, cs)

    assert exons == [(1001, 1103)]
    assert queries == [(0, 97)]
    assert mismatches == 0
    assert gaps == 1
    assert proper_splice is True


@patch("gapmm2.align.execute")
def test_splice_aligner_minimap2_missing_tags(mock_execute):
    """Test splice_aligner_minimap2 skips alignments with missing optional tags."""
    # Mock minimap2 output with missing ts tag (which caused the original bug)
    mock_output = [
        "query1\t100\t0\t100\tref1\t1000\t500\t600\t100\t100\t60\t255\ttp:A:P\tNM:i:0\tcs:Z:=100"
    ]
    mock_execute.return_value = mock_output

    # This should not raise UnboundLocalError anymore
    results = list(splice_aligner_minimap2("ref.fa", "query.fa", refine=False))

    # Should get only stats (no alignment results since ts tag is missing)
    assert len(results) == 1

    # Check stats - should track missing ts tag and skip alignment
    stats = results[0]
    assert stats["n"] == 1  # One alignment processed
    assert stats["low-mapq"] == 0
    assert stats["unaligned"]["missing-ts"] == 1  # ts tag was missing, alignment skipped
    assert stats["unaligned"]["missing-tp"] == 0  # tp tag was present
    assert stats["unaligned"]["missing-nm"] == 0  # nm tag was present
    assert stats["unaligned"]["missing-cs"] == 0  # cs tag was present


@patch("gapmm2.align.execute")
def test_splice_aligner_minimap2_missing_cs_tag(mock_execute):
    """Test splice_aligner_minimap2 skips alignments missing required cs tag."""
    # Mock minimap2 output with no cs tag (should be skipped)
    mock_output = ["query1\t100\t0\t100\tref1\t1000\t500\t600\t100\t100\t60\t255\ttp:A:P\tNM:i:0"]
    mock_execute.return_value = mock_output

    # This should not raise UnboundLocalError anymore
    results = list(splice_aligner_minimap2("ref.fa", "query.fa", refine=False))

    # Should get only stats (no alignment results since cs tag is missing)
    assert len(results) == 1

    # Check stats - should track missing cs tag and skip alignment
    stats = results[0]
    assert stats["n"] == 1  # One alignment processed
    assert stats["unaligned"]["missing-cs"] == 1  # cs tag was missing, alignment skipped
    assert stats["unaligned"]["invalid-paf"] == 0  # PAF format was valid


@patch("gapmm2.align.execute")
def test_splice_aligner_minimap2_all_optional_tags_missing(mock_execute):
    """Test splice_aligner_minimap2 skips alignments with missing optional tags."""
    # Mock minimap2 output with only cs tag present (missing tp, ts, nm)
    mock_output = ["query1\t100\t0\t100\tref1\t1000\t500\t600\t100\t100\t60\t255\tcs:Z:=100"]
    mock_execute.return_value = mock_output

    # This should not raise UnboundLocalError anymore
    results = list(splice_aligner_minimap2("ref.fa", "query.fa", refine=False))

    # Should get only stats (no alignment results since tp tag is missing and checked first)
    assert len(results) == 1

    # Check stats - should track missing tp tag (first one checked)
    stats = results[0]
    assert stats["n"] == 1
    assert stats["low-mapq"] == 0
    assert stats["unaligned"]["missing-tp"] == 1  # tp tag was missing, alignment skipped
    assert stats["unaligned"]["missing-ts"] == 0  # ts check not reached
    assert stats["unaligned"]["missing-nm"] == 0  # nm check not reached
    assert stats["unaligned"]["missing-cs"] == 0  # cs tag was present


@patch("gapmm2.align.execute")
def test_splice_aligner_minimap2_invalid_paf(mock_execute):
    """Test splice_aligner_minimap2 handles invalid PAF lines gracefully."""
    # Mock minimap2 output with invalid PAF (too few columns)
    mock_output = [
        "query1\t100\t0\t100\tref1",  # Only 5 columns, should have at least 12
        "",  # Empty line
        "query2\t200\t0\t200\tref1\t2000\t1500\t1700\t200\t200\t60\t255\ttp:A:P\tts:A:+\tNM:i:0\tcs:Z:=200",  # Valid line with all tags
    ]
    mock_execute.return_value = mock_output

    # This should not raise errors
    results = list(splice_aligner_minimap2("ref.fa", "query.fa", refine=False))

    # Should get one alignment result plus stats
    assert len(results) == 2

    # Check that valid alignment was processed
    paf_result = results[0]
    assert len(paf_result) == 16  # 12 standard PAF fields + 4 extra fields
    assert paf_result[12] == "tp:A:P"  # tp tag present
    assert paf_result[13] == "ts:A:+"  # ts tag present
    assert paf_result[14] == "NM:i:0"  # nm tag present
    assert paf_result[15] == "cs:Z:=200"  # cs tag present

    # Check stats
    stats = results[1]
    assert stats["n"] == 2  # Two lines processed (empty line not counted)
    assert stats["unaligned"]["invalid-paf"] == 1  # One invalid PAF line
    assert stats["unaligned"]["missing-cs"] == 0  # Valid line had cs tag


@patch("gapmm2.align.execute")
def test_splice_aligner_minimap2_successful_alignment(mock_execute):
    """Test splice_aligner_minimap2 with successful alignment containing all tags."""
    # Mock minimap2 output with all required tags present
    mock_output = [
        "query1\t100\t0\t100\tref1\t1000\t500\t600\t100\t100\t60\t255\ttp:A:P\tts:A:+\tNM:i:0\tcs:Z:=100"
    ]
    mock_execute.return_value = mock_output

    # This should work without issues
    results = list(splice_aligner_minimap2("ref.fa", "query.fa", refine=False))

    # Should get one alignment result plus stats
    assert len(results) == 2

    # Check that the PAF format is correct
    paf_result = results[0]
    assert len(paf_result) == 16  # 12 standard PAF fields + 4 extra fields
    assert paf_result[12] == "tp:A:P"  # tp tag present
    assert paf_result[13] == "ts:A:+"  # ts tag present
    assert paf_result[14] == "NM:i:0"  # nm tag present
    assert paf_result[15] == "cs:Z:=100"  # cs tag present

    # Check stats - should have no unaligned reads
    stats = results[1]
    assert stats["n"] == 1
    assert stats["low-mapq"] == 0
    assert stats["unaligned"]["missing-tp"] == 0
    assert stats["unaligned"]["missing-ts"] == 0
    assert stats["unaligned"]["missing-nm"] == 0
    assert stats["unaligned"]["missing-cs"] == 0
    assert stats["unaligned"]["invalid-paf"] == 0


@patch("gapmm2.align.mp.Aligner")
@patch("gapmm2.align.mp.fastx_read")
def test_splice_aligner_nested_stats(mock_fastx_read, mock_aligner):
    """Test splice_aligner has the same nested statistics structure."""
    # Mock the mappy components
    mock_aligner_instance = MagicMock()
    mock_aligner.return_value = mock_aligner_instance

    # Mock alignment hit with all required attributes
    mock_hit = MagicMock()
    mock_hit.mapq = 60
    mock_hit.strand = 1
    mock_hit.cs = "=100"
    mock_hit.NM = 0
    mock_hit.is_primary = 1
    mock_hit.trans_strand = 1
    mock_hit.ctg = "ref1"
    mock_hit.ctg_len = 1000
    mock_hit.r_st = 500
    mock_hit.r_en = 600
    mock_hit.q_st = 0
    mock_hit.q_en = 100
    mock_hit.mlen = 100
    mock_hit.blen = 100

    mock_aligner_instance.map.return_value = [mock_hit]
    mock_fastx_read.return_value = [("query1", "A" * 100, None)]

    # Run the function
    results = list(splice_aligner("ref.fa", "query.fa", refine=False))

    # Should get one alignment result plus stats
    assert len(results) == 2

    # Check that the PAF format is correct
    paf_result = results[0]
    assert len(paf_result) == 16  # 12 standard PAF fields + 4 extra fields

    # Check stats - should have nested unaligned structure
    stats = results[1]
    assert stats["n"] == 1
    assert stats["low-mapq"] == 0
    assert "unaligned" in stats
    assert stats["unaligned"]["missing-tp"] == 0
    assert stats["unaligned"]["missing-ts"] == 0
    assert stats["unaligned"]["missing-nm"] == 0
    assert stats["unaligned"]["missing-cs"] == 0
    assert stats["unaligned"]["invalid-paf"] == 0
