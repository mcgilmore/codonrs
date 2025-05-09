# codonrs

`codonrs` is a small crate for rapidly calculating relative synonymous codon usage (RSCU)
values for coding DNA sequences, for analyses of codon usage bias. The crate can be used
as a command-line utility with the `codonrs` command, or used in other crates via the analysis mod.

### Command-line usage

#### Required arguments

**Input**: `-i`/`--input`: A multi-fasta file with the sequences to be analysed as individual
fasta entries. Sequences whose length is not a multiple of three will be ignored

**Output**: `-o`/`--output`: Prefix for the output files. Three files are output by default:
  - `prefix`_codon.csv: raw codon counts for each CDS.
  - `prefix`_amino_acids.csv: amino acid counts for each CDS, determined from the chosen translation table.
  - `prefix`_rscu.csv: calculated RSCU values for each CDS.

#### Optional arguments

**Translation table**: `-t`/`--table`: Integer representing NCBI translation table to be used for codon
counts and RSCU calculation. See tables [here](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes).
Defaults to 1: the standard code.



License: GPL-3.0
