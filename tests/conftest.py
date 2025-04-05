#!/usr/bin/env python3

import os
import pytest
import tempfile


@pytest.fixture
def temp_fasta_files():
    """Create temporary FASTA files for testing."""
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a reference genome FASTA file
        ref_path = os.path.join(tmpdir, "reference.fa")
        with open(ref_path, "w") as f:
            f.write(">scaffold_1\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
        
        # Create a query transcript FASTA file
        query_path = os.path.join(tmpdir, "query.fa")
        with open(query_path, "w") as f:
            f.write(">transcript_1\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
        
        yield ref_path, query_path
