[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/gapmm2.svg)](https://github.com/nextgenusfs/gapmm2/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/gapmm2)
[![codecov](https://codecov.io/gh/nextgenusfs/gapmm2/graph/badge.svg?token=THI1LXFFNX)](https://codecov.io/gh/nextgenusfs/gapmm2)

# gapmm2: gapped alignment using minimap2

This tool is a wrapper for minimap2 to run spliced/gapped alignment, ie aligning transcripts to a genome.   You are probably saying, yes minimap2 runs this with `-x splice --cs` option (you are correct).  However, there are instances where the terminal exons from stock minimap2 alignments are missing. This tool detects those alignments that have unaligned terminal eons and uses `edlib` to find the terminal exon positions. The tool then updates the PAF output file with the updated information.

#### Rationale

We can pull out a gene model in GFF3 format that has a short 5' terminal exon:

```
scaffold_9	funannotate	gene	408904	409621	.	-	.	ID=OPO1_006919;
scaffold_9	funannotate	mRNA	408904	409621	.	-	.	ID=OPO1_006919-T1;Parent=OPO1_006919;product=hypothetical protein;
scaffold_9	funannotate	exon	409609	409621	.	-	.	ID=OPO1_006919-T1.exon1;Parent=OPO1_006919-T1;
scaffold_9	funannotate	exon	409320	409554	.	-	.	ID=OPO1_006919-T1.exon2;Parent=OPO1_006919-T1;
scaffold_9	funannotate	exon	409090	409255	.	-	.	ID=OPO1_006919-T1.exon3;Parent=OPO1_006919-T1;
scaffold_9	funannotate	exon	408904	409032	.	-	.	ID=OPO1_006919-T1.exon4;Parent=OPO1_006919-T1;
scaffold_9	funannotate	CDS	409609	409621	.	-	0	ID=OPO1_006919-T1.cds;Parent=OPO1_006919-T1;
scaffold_9	funannotate	CDS	409320	409554	.	-	2	ID=OPO1_006919-T1.cds;Parent=OPO1_006919-T1;
scaffold_9	funannotate	CDS	409090	409255	.	-	1	ID=OPO1_006919-T1.cds;Parent=OPO1_006919-T1;
scaffold_9	funannotate	CDS	408904	409032	.	-	0	ID=OPO1_006919-T1.cds;Parent=OPO1_006919-T1;
```

If we then map this transcript against the genome, we get the following PAF alignment.

```
$ minimap2 -x splice --cs genome.fasta cds-transcripts.fa | grep 'OPO1_006919'
OPO1_006919-T1	543	13	543	-	scaffold_9	658044	408903	409554	530	530	60	NM:i:0	ms:i:530	AS:i:466	nn:i:0	ts:A:+	tp:A:P	cm:i:167	s1:i:510	s2:i:0	de:f:0	rl:i:0	cs:Z::129~ct57ac:166~ct64ac:235
```

The `--cs` flag in minimap2 can be used to parse the coordinates (below) and you can see we are missing the 5' exon.

```
>>> cs2coords(408903, 13, 543, '-', ':129~ct57ac:166~ct64ac:235')
([(409320, 409554), (409090, 409255), (408904, 409032)],
```

So if we run this same alignment with `gapmm2` we are able to properly align the 5' terminal exon.

```
$ gapmm2 genome.fa cds-transcripts.fa | grep 'OPO1_006919'
OPO1_006919-T1	543	0	543	-	scaffold_9	658044	408903	409621	543	543	60	tp:A:P	ts:A:+	NM:i:0	cs:Z::129~ct57ac:166~ct64ac:235~ct54ac:13
```

```
>>> cs2coords(408903, 0, 543, '-', ':129~ct57ac:166~ct64ac:235~ct54ac:13')
([(409609, 409621), (409320, 409554), (409090, 409255), (408904, 409032)]
```



#### Usage:

`gapmm2` can be run as a command line script:

```
$ gapmm2
usage: gapmm2 [-o] [-f] [-t] [-m] [-i] [-d] [-h] [--version] reference query

gapmm2: gapped alignment with minimap2. Performs minimap2/mappy alignment with splice options and refines terminal alignments with edlib.

Positional arguments:
  reference           reference genome (FASTA)
  query               transcipts in FASTA or FASTQ

Optional arguments:
  -o , --out          output in PAF format (default: stdout)
  -f , --out-format   output format [paf,gff3] (default: paf)
  -t , --threads      number of threads to use with minimap2 (default: 3)
  -m , --min-mapq     minimum map quality value (default: 1)
  -i , --max-intron   max intron length, controls terminal search space (default: 500)
  -d, --debug         write some debug info to stderr (default: False)

Help:
  -h, --help          Show this help message and exit
  --version           Show program's version number and exit
```



### Python API

It can also be run as a python module. The module provides several functions for working with spliced alignments:

#### `splice_aligner` function

The main function for aligning transcripts to a genome. It returns a list of lists containing PAF formatted data for each alignment and a dictionary of simple stats.

```python
>>> from gapmm2.align import splice_aligner
>>> results, stats = splice_aligner('genome.fa', 'transcripts.fa')
>>> stats
{'n': 6926, 'low-mapq': 0, 'refine-left': 409, 'refine-right': 63}
>>> len(results)
6926
>>> results[0]
['OPO1_000001-T1', 2184, 0, 2184, '+', 'scaffold_1', 1803704, 887, 3127, 2184, 2184, 60, 'tp:A:P', 'ts:A:+', 'NM:i:0', 'cs:Z::958~gt56ag:1226']
```

#### `cs2coords` function

This function parses the CIGAR string (cs) from minimap2 and converts it to genomic coordinates, identifying exons, introns, and other alignment features.

```python
>>> from gapmm2.align import cs2coords
>>> cs2coords(408903, 0, 543, '-', ':129~ct57ac:166~ct64ac:235~ct54ac:13')
([(409609, 409621), (409320, 409554), (409090, 409255), (408904, 409032)], [(0, 13), (13, 248), (248, 414), (414, 543)], 0, 0, True)
```

#### `paf2gff3` function

This function converts PAF format alignments to GFF3 format.

```python
>>> from gapmm2.align import splice_aligner, paf2gff3
>>> results, stats = splice_aligner('genome.fa', 'transcripts.fa')
>>> paf2gff3(results, output='output.gff3')
```



### Installation

You can install gapmm2 using pip:

```bash
pip install gapmm2
```

Or you can install the latest development version directly from GitHub:

```bash
pip install git+https://github.com/nextgenusfs/gapmm2.git
```

You can also install from conda:

```bash
conda install -c bioconda gapmm2
```

### Dependencies

Gapmm2 requires the following Python packages:

- mappy (Python bindings for minimap2)
- edlib (for sequence alignment)
- natsort (for natural sorting)

These dependencies will be automatically installed when you install gapmm2 using pip or conda.

### Development

#### Testing

Gapmm2 includes a test suite that can be run using pytest. To run the tests, first install pytest:

```bash
pip install pytest pytest-cov
```

Then run the tests from the root directory of the repository:

```bash
python -m pytest tests/ --cov=gapmm2
```

#### Code Formatting

This project uses [pre-commit](https://pre-commit.com/) to ensure code quality and consistency. The pre-commit hooks run Black (code formatter), isort (import sorter), and flake8 (linter).

To set up pre-commit:

1. Install pre-commit:

```bash
pip install pre-commit
```

2. Install the git hooks:

```bash
pre-commit install
```

3. (Optional) Run against all files:

```bash
pre-commit run --all-files
```

After installation, the pre-commit hooks will run automatically on each commit to ensure your code follows the project's style guidelines.
