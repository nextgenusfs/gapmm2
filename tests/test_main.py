#!/usr/bin/env python3

import sys
import pytest
from io import StringIO
from unittest.mock import patch
from gapmm2.__main__ import main


def test_main_help():
    """Test main function with --help argument."""
    # Mock sys.argv and capture the output
    with (
        patch.object(sys, "argv", ["gapmm2", "--help"]),
        patch("sys.stdout", new_callable=StringIO) as mock_stdout,
        pytest.raises(SystemExit) as excinfo,
    ):
        main()

    # Check that the exit code is 0 (success)
    assert excinfo.value.code == 0

    # Check that the help text contains expected information
    output = mock_stdout.getvalue()
    assert "gapmm2: gapped alignment with minimap2" in output
    assert "reference" in output
    assert "query" in output


def test_main_version():
    """Test main function with --version argument."""
    # Mock sys.argv and capture the output
    with (
        patch.object(sys, "argv", ["gapmm2", "--version"]),
        patch("sys.stdout", new_callable=StringIO) as mock_stdout,
        pytest.raises(SystemExit) as excinfo,
    ):
        main()

    # Check that the exit code is 0 (success)
    assert excinfo.value.code == 0

    # Check that the version is printed
    output = mock_stdout.getvalue()
    assert "gapmm2 v" in output
