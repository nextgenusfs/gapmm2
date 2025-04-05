#!/usr/bin/env python3

import pytest
from gapmm2.utils import check_inputs, zopen


def test_check_inputs_valid_files(tmp_path):
    """Test check_inputs with valid files."""
    # Create temporary files
    file1 = tmp_path / "file1.txt"
    file2 = tmp_path / "file2.txt"
    file1.write_text("test content 1")
    file2.write_text("test content 2")

    # This should not raise an exception
    check_inputs([str(file1), str(file2)])


def test_check_inputs_invalid_files(tmp_path):
    """Test check_inputs with invalid files."""
    # Create one valid file
    file1 = tmp_path / "file1.txt"
    file1.write_text("test content 1")

    # Non-existent file
    file2 = tmp_path / "nonexistent.txt"

    # This should raise a FileNotFoundError exception
    with pytest.raises(FileNotFoundError):
        check_inputs([str(file1), str(file2)])


def test_zopen_read_text(tmp_path):
    """Test zopen for reading text files."""
    # Create a text file
    test_file = tmp_path / "test.txt"
    test_file.write_text("test content")

    # Open and read the file
    with zopen(str(test_file), mode="r") as f:
        content = f.read()

    assert content == "test content"


def test_zopen_write_text(tmp_path):
    """Test zopen for writing text files."""
    # Define a file path
    test_file = tmp_path / "output.txt"

    # Write to the file
    with zopen(str(test_file), mode="w") as f:
        f.write("test output")

    # Verify the content
    assert test_file.read_text() == "test output"
