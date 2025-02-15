use clap::Parser;
use std::collections::{HashMap, HashSet};
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
fn write_rscu_to_csv(
    filename: &str,
    rscu_data: &Vec<(String, HashMap<String, f64>)>,
) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    // Collect all unique codons across sequences to create CSV headers
    let mut codon_set = HashSet::new();
    for (_, codon_map) in rscu_data {
        for codon in codon_map.keys() {
            codon_set.insert(codon.clone());
        }
    }

    let mut codons: Vec<String> = codon_set.into_iter().collect();
    codons.sort(); // Ensure consistent order

    // Write the CSV header
    write!(file, "Sequence")?;
    for codon in &codons {
        write!(file, ",{}", codon)?;
    }
    writeln!(file)?;

    // Write RSCU values for each sequence
    for (seq_name, codon_map) in rscu_data {
        write!(file, "{}", seq_name)?;
        for codon in &codons {
            let rscu_value = codon_map.get(codon).unwrap_or(&0.0);
            write!(file, ",{:.6}", rscu_value)?;
        }
        writeln!(file)?;
    }

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

/// Writes codon and amino acid counts to a CSV file.
fn write_counts_to_csv(
    filename_prefix: &str,
    codon_data: &Vec<(String, HashMap<String, usize>)>,
    amino_acid_data: &Vec<(String, HashMap<String, usize>)>,
) -> std::io::Result<()> {
    let codon_filename = format!("{}_codon.csv", filename_prefix);
    let amino_filename = format!("{}_amino_acids.csv", filename_prefix);

    // Create files
    let mut codon_file = File::create(&codon_filename)?;
    let mut amino_file = File::create(&amino_filename)?;

    // Collect all unique codons and amino acids across sequences
    let mut codon_set = HashSet::new();
    let mut amino_acid_set = HashSet::new();

    for (_, codon_map) in codon_data {
        for codon in codon_map.keys() {
            codon_set.insert(codon.clone());
        }
    }

    for (_, amino_map) in amino_acid_data {
        for amino in amino_map.keys() {
            amino_acid_set.insert(amino.clone());
        }
    }

    let mut codons: Vec<String> = codon_set.into_iter().collect();
    let mut amino_acids: Vec<String> = amino_acid_set.into_iter().collect();
    codons.sort(); // Ensure consistent order
    amino_acids.sort();

    // Write the CSV header for codons
    write!(codon_file, "Sequence")?;
    for codon in &codons {
        write!(codon_file, ",{}", codon)?;
    }
    writeln!(codon_file)?;

    // Write codon counts for each sequence
    for (seq_name, codon_map) in codon_data {
        write!(codon_file, "{}", seq_name)?;
        for codon in &codons {
            let count = codon_map.get(codon).unwrap_or(&0);
            write!(codon_file, ",{}", count)?;
        }
        writeln!(codon_file)?;
    }

    // Write the CSV header for amino acids
    write!(amino_file, "Sequence")?;
    for amino in &amino_acids {
        write!(amino_file, ",{}", amino)?;
    }
    writeln!(amino_file)?;

    // Write amino acid counts for each sequence
    for (seq_name, amino_map) in amino_acid_data {
        write!(amino_file, "{}", seq_name)?;
        for amino in &amino_acids {
            let count = amino_map.get(amino).unwrap_or(&0);
            write!(amino_file, ",{}", count)?;
        }
        writeln!(amino_file)?;
    }

    Ok(())
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
    let args = Cli::parse();
    let mut rscu_results = Vec::new();
    let mut codon_counts_list = Vec::new();
    let mut amino_acid_counts_list = Vec::new();

    match read_sequences_from_fasta(&args.input_file) {
        Ok(sequences) => {
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

                // Store data for CSV output
                rscu_results.push((seq_name.clone(), rscu_values));
                codon_counts_list.push((seq_name.clone(), codon_counts));
                amino_acid_counts_list.push((seq_name.clone(), amino_acid_counts));
            }

            // Write RSCU values to a single CSV file
            let rscu_filename = format!("{}_rscu.csv", args.output_file);
            write_rscu_to_csv(&rscu_filename, &rscu_results).unwrap();

            // Write Codon & Amino Acid Counts to CSV
            write_counts_to_csv(
                &args.output_file,
                &codon_counts_list,
                &amino_acid_counts_list,
            )
            .unwrap();

            println!(
                "Results saved to:\n - {}_counts.csv\n - {}_amino_acids.csv\n - {}_rscu.csv",
                &args.output_file, &args.output_file, &args.output_file,
            );
        }
        Err(e) => eprintln!("Error reading file: {}", e),
    }
}
