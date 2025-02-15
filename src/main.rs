use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// Command-line arguments using Clap
#[derive(Parser)]
#[command(name = "CodonRS")]
#[command(version = "0.1.0")]
#[command(about = "Analyze codon usage bias in DNA sequences", long_about = None)]
struct Cli {
    /// Input DNA sequence file (FASTA or plain text)
    #[arg(short = 'i', long = "input")]
    input_file: String,

    /// Output file for results
    #[arg(short = 'o', long = "output")]
    output_file: String,

    /// NCBI translation table ID (default: 1)
    #[arg(short = 't', long = "table", default_value_t = 1)]
    translation_table: u8,
}

// Function to load NCBI translation tables
fn get_ncbi_translation_table(table_id: u8) -> HashMap<&'static str, &'static str> {
    match table_id {
        1 => HashMap::from([
            ("TTT", "F"),
            ("TTC", "F"),
            ("TTA", "L"),
            ("TTG", "L"),
            ("TCT", "S"),
            ("TCC", "S"),
            ("TCA", "S"),
            ("TCG", "S"),
            ("TAT", "Y"),
            ("TAC", "Y"),
            ("TAA", "*"),
            ("TAG", "*"),
            ("TGT", "C"),
            ("TGC", "C"),
            ("TGA", "*"),
            ("TGG", "W"),
            ("CTT", "L"),
            ("CTC", "L"),
            ("CTA", "L"),
            ("CTG", "L"),
            ("CCT", "P"),
            ("CCC", "P"),
            ("CCA", "P"),
            ("CCG", "P"),
            ("CAT", "H"),
            ("CAC", "H"),
            ("CAA", "Q"),
            ("CAG", "Q"),
            ("CGT", "R"),
            ("CGC", "R"),
            ("CGA", "R"),
            ("CGG", "R"),
            ("ATT", "I"),
            ("ATC", "I"),
            ("ATA", "I"),
            ("ATG", "M"),
            ("ACT", "T"),
            ("ACC", "T"),
            ("ACA", "T"),
            ("ACG", "T"),
            ("AAT", "N"),
            ("AAC", "N"),
            ("AAA", "K"),
            ("AAG", "K"),
            ("AGT", "S"),
            ("AGC", "S"),
            ("AGA", "R"),
            ("AGG", "R"),
            ("GTT", "V"),
            ("GTC", "V"),
            ("GTA", "V"),
            ("GTG", "V"),
            ("GCT", "A"),
            ("GCC", "A"),
            ("GCA", "A"),
            ("GCG", "A"),
            ("GAT", "D"),
            ("GAC", "D"),
            ("GAA", "E"),
            ("GAG", "E"),
            ("GGT", "G"),
            ("GGC", "G"),
            ("GGA", "G"),
            ("GGG", "G"),
        ]),
        _ => panic!("Translation table {} not implemented!", table_id),
    }
}

/// Parse codon content in sequence
fn parse_codons(sequence: &str) -> HashMap<String, usize> {
    let mut codon_counts = HashMap::new();

    // Ignore sequences that are not a multiple of 3
    if sequence.len() % 3 != 0 {
        return codon_counts; // Return an empty HashMap instead of None
    }

    for codon in sequence.as_bytes().chunks(3) {
        let codon_str = String::from_utf8_lossy(codon).to_uppercase();
        *codon_counts.entry(codon_str).or_insert(0) += 1;
    }

    codon_counts
}

fn compute_rscu(
    codon_counts: &HashMap<String, usize>,
    translation_table: u8,
) -> HashMap<String, f64> {
    let codon_table = get_ncbi_translation_table(translation_table);
    let mut rscu_values = HashMap::new();
    let mut amino_acid_totals: HashMap<&str, usize> = HashMap::new();
    let mut synonymous_codons: HashMap<&str, Vec<&str>> = HashMap::new();

    // Group codons by their amino acid and count occurrences
    for (codon, amino_acid) in &codon_table {
        synonymous_codons
            .entry(amino_acid)
            .or_insert(Vec::new())
            .push(codon);
        *amino_acid_totals.entry(amino_acid).or_insert(0) +=
            codon_counts.get(*codon).copied().unwrap_or(0);
    }

    // Compute RSCU values
    for (amino_acid, codons) in &synonymous_codons {
        let total_codon_count = amino_acid_totals.get(amino_acid).copied().unwrap_or(0) as f64;
        let num_codons = codons.len() as f64;

        for codon in codons {
            let observed = *codon_counts.get(*codon).unwrap_or(&0) as f64;
            let expected = total_codon_count / num_codons;
            let rscu = if expected > 0.0 {
                observed / expected
            } else {
                0.0
            };

            rscu_values.insert((*codon).to_string(), rscu);
        }
    }

    rscu_values
}

/// Writes RSCU values to a CSV file.
fn write_rscu_to_csv(filename: &str, rscu_data: &HashMap<String, f64>) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    writeln!(file, "Codon,RSCU")?; // CSV Header

    for (codon, rscu) in rscu_data {
        writeln!(file, "{},{}", codon, rscu)?;
    }

    println!("RSCU values saved to {}", filename);
    Ok(())
}

/// Translate DNA sequence into amino acids
fn translate_sequence(sequence: &str, table_id: u8) -> HashMap<String, usize> {
    let codon_table = get_ncbi_translation_table(table_id);
    let mut amino_acid_counts = HashMap::new();

    for codon in sequence.as_bytes().chunks(3) {
        if codon.len() == 3 {
            let codon_str = String::from_utf8_lossy(codon).to_uppercase();
            if let Some(amino_acid) = codon_table.get(codon_str.as_str()) {
                *amino_acid_counts.entry(amino_acid.to_string()).or_insert(0) += 1;
            }
        }
    }

    amino_acid_counts
}

/// Read sequences from a multi-FASTA file
fn read_sequences_from_fasta(filename: &str) -> io::Result<Vec<(String, String)>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut current_seq_name = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_seq_name.is_empty() {
                sequences.push((current_seq_name.clone(), current_sequence.clone()));
                current_sequence.clear();
            }
            current_seq_name = line[1..].to_string();
        } else {
            current_sequence.push_str(&line);
        }
    }
    if !current_seq_name.is_empty() {
        sequences.push((current_seq_name, current_sequence));
    }

    Ok(sequences)
}

fn main() {
    let args = Cli::parse(); // Parse command-line arguments

    match read_sequences_from_fasta(&args.input_file) {
        Ok(sequences) => {
            let mut output_file =
                File::create(&args.output_file).expect("Failed to create output file");

            for (seq_name, sequence) in sequences {
                let codon_counts = parse_codons(&sequence);

                if codon_counts.is_empty() {
                    eprintln!(
                        "Warning: Skipping sequence '{}' (not a multiple of 3)",
                        seq_name
                    );
                    continue;
                }

                let amino_acid_counts = translate_sequence(&sequence, args.translation_table);
                let rscu_values = compute_rscu(&codon_counts, args.translation_table);

                writeln!(output_file, ">{}", seq_name).unwrap();
                writeln!(output_file, "Codon Composition:").unwrap();
                for (codon, count) in &codon_counts {
                    writeln!(output_file, "{}: {}", codon, count).unwrap();
                }

                writeln!(output_file, "\nAmino Acid Composition:").unwrap();
                for (amino, count) in &amino_acid_counts {
                    writeln!(output_file, "{}: {}", amino, count).unwrap();
                }
                writeln!(output_file, "\n").unwrap();

                // Write RSCU values to a CSV file
                let rscu_filename = format!("{}_rscu.csv", seq_name);
                write_rscu_to_csv(&rscu_filename, &rscu_values).unwrap();
            }

            println!("Results saved to {}", args.output_file);
        }
        Err(e) => eprintln!("Error reading file: {}", e),
    }
}
