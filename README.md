# codonrs

Calculate relative synonymous codon usage (RSCU) for coding DNA sequences (CDS).

## Usage

`codonrs [Optional: -t  <translation table>] -i <input.fasta> -o <prefix for output CSV files>`

- `<translation table>` - (Optional) an integer representing the translation table to use, following [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes) numbering. Defaults to 1.
- `input.fasta` - a multi fasta file containing DNA sequences to be analysed.
