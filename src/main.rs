use clap::Parser;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

/// Load embedded JSON file at compile-time
static CODE_FILE: &str = include_str!("genetic_code.json");

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

    /// Flag to compute and output Z-score
    #[arg(short = 'z', long = "zscore", default_value_t = false)]
    compute_zscore: bool,
}

#[derive(Debug, Clone)]
struct GeneticCode {
    id: String,
    name: String,
    codon_map: HashMap<String, String>,
}

/// A genetic code object, containing an id, name and translation table (codon_map)
impl GeneticCode {
    fn new() -> Self {
        GeneticCode {
            id: String::new(),
            name: String::new(),
            codon_map: HashMap::new(),
        }
    }
}

/// Define the structure for the JSON format
#[derive(Debug, Serialize, Deserialize)]
struct GeneticCodeJSON {
    name: Vec<String>,
    id: u8,
    ncbieaa: String,
    sncbieaa: String,
    base_mappings: HashMap<String, String>,
}

impl GeneticCodeJSON {
    /// Convert the genetic code into a codon-to-amino-acid mapping
    fn to_codon_map(&self) -> HashMap<String, String> {
        let base1 = self
            .base_mappings
            .get("Base1")
            .cloned()
            .unwrap_or_else(|| String::new());
        let base2 = self
            .base_mappings
            .get("Base2")
            .cloned()
            .unwrap_or_else(|| String::new());
        let base3 = self
            .base_mappings
            .get("Base3")
            .cloned()
            .unwrap_or_else(|| String::new());

        let mut codon_map = HashMap::new();

        for (i, ((b1, b2), b3)) in base1
            .chars()
            .zip(base2.chars())
            .zip(base3.chars())
            .enumerate()
        {
            let codon = format!("{}{}{}", b1, b2, b3);
            let amino_acid = self.ncbieaa.chars().nth(i).unwrap_or('*').to_string(); // Convert char to String
            codon_map.insert(codon, amino_acid);
        }

        codon_map
    }
}

/// Load the provided genetic_code.json data file
fn load_embedded_genetic_codes() -> Result<Vec<GeneticCodeJSON>, Box<dyn Error>> {
    let genetic_codes: Vec<GeneticCodeJSON> = serde_json::from_str(CODE_FILE)?;
    Ok(genetic_codes)
}

/// Get a genetic code from a GeneticCode object by a provided ID
fn get_genetic_code_by_id<'a>(
    codes: &'a [GeneticCodeJSON],
    id: &u8,
) -> Option<&'a GeneticCodeJSON> {
    codes.iter().find(|code| code.id == *id)
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

fn compute_rscu(codon_counts: &HashMap<String, usize>, code: &GeneticCode) -> HashMap<String, f64> {
    let codon_table = &code.codon_map;
    let mut rscu_values = HashMap::new();
    let mut amino_acid_totals: HashMap<&str, usize> = HashMap::new();
    let mut synonymous_codons: HashMap<&str, Vec<&str>> = HashMap::new();

    // Group codons by their amino acid and count occurrences
    for (codon, amino_acid) in codon_table {
        synonymous_codons
            .entry(amino_acid)
            .or_insert(Vec::new())
            .push(codon);
        *amino_acid_totals.entry(amino_acid).or_insert(0) +=
            codon_counts.get(codon).copied().unwrap_or(0);
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

/// Compute the mean RSCU value for each codon
fn compute_mean_rscu(rscu_results: &Vec<(String, HashMap<String, f64>)>) -> HashMap<String, f64> {
    let mut total_rscu: HashMap<String, f64> = HashMap::new();
    let mut gene_count = rscu_results.len() as f64;

    for (_, rscu_map) in rscu_results {
        for (codon, value) in rscu_map {
            *total_rscu.entry(codon.clone()).or_insert(0.0) += value;
        }
    }

    let mut mean_rscu = HashMap::new();
    for (codon, total_value) in total_rscu {
        mean_rscu.insert(codon, total_value / gene_count);
    }

    mean_rscu
}

/// Compute standard deviation of RSCU values
fn compute_std_rscu(
    rscu_results: &Vec<(String, HashMap<String, f64>)>,
    mean_rscu: &HashMap<String, f64>,
) -> HashMap<String, f64> {
    let mut variance: HashMap<String, f64> = HashMap::new();
    let gene_count = rscu_results.len() as f64;

    for (_, rscu_map) in rscu_results {
        for (codon, value) in rscu_map {
            let mean = mean_rscu.get(codon).unwrap_or(&0.0);
            let diff = value - mean;
            *variance.entry(codon.clone()).or_insert(0.0) += diff * diff;
        }
    }

    let mut std_dev = HashMap::new();
    for (codon, var) in variance {
        std_dev.insert(codon, (var / gene_count).sqrt());
    }

    std_dev
}

/// Compute Z-score from RSCU values
fn compute_rscu_z_scores(
    rscu_results: &Vec<(String, HashMap<String, f64>)>,
    mean_rscu: &HashMap<String, f64>,
    std_rscu: &HashMap<String, f64>,
) -> Vec<(String, HashMap<String, f64>)> {
    let mut z_scores = Vec::new();

    for (gene, rscu_map) in rscu_results {
        let mut gene_z_scores = HashMap::new();
        for (codon, value) in rscu_map {
            let mean = mean_rscu.get(codon).unwrap_or(&0.0);
            let std_dev = std_rscu.get(codon).unwrap_or(&1.0); // Avoid division by zero
            let z_score = (value - mean) / std_dev;
            gene_z_scores.insert(codon.clone(), z_score);
        }
        z_scores.push((gene.clone(), gene_z_scores));
    }

    z_scores
}

/// Write Z-scores to a CSV file
fn write_z_scores_to_csv(
    filename: &str,
    z_scores: &Vec<(String, HashMap<String, f64>)>,
) -> io::Result<()> {
    let mut file = File::create(filename)?;

    let mut codon_set = HashSet::new();
    for (_, codon_map) in z_scores {
        for codon in codon_map.keys() {
            codon_set.insert(codon.clone());
        }
    }

    let mut codons: Vec<String> = codon_set.into_iter().collect();
    codons.sort();

    // Write header
    write!(file, "Gene")?;
    for codon in &codons {
        write!(file, ",{}", codon)?;
    }
    writeln!(file)?;

    // Write Z-scores for each gene
    for (gene, codon_map) in z_scores {
        write!(file, "{}", gene)?;
        for codon in &codons {
            let z_score = codon_map.get(codon).unwrap_or(&0.0);
            write!(file, ",{:.6}", z_score)?;
        }
        writeln!(file)?;
    }

    Ok(())
}

/// Translate DNA sequence into amino acids
fn translate_sequence(sequence: &str, code: &GeneticCode) -> HashMap<String, usize> {
    let codon_table = &code.codon_map;
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

    let mut code = GeneticCode::new();

    match load_embedded_genetic_codes() {
        Ok(genetic_codes) => {
            if let Some(selected_code) =
                get_genetic_code_by_id(&genetic_codes, &args.translation_table)
            {
                code.id = selected_code.id.to_string();
                code.name = selected_code
                    .name
                    .first()
                    .cloned()
                    .unwrap_or_else(|| String::new());
                code.codon_map = selected_code.to_codon_map();
                println!("Using genetic code {}: {}", code.id, code.name);
            } else {
                eprintln!(
                    "Error: Genetic Code ID {} not found!",
                    &args.translation_table
                );
            }
        }
        Err(e) => eprintln!("Failed to load genetic codes: {}", e),
    }

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

                let amino_acid_counts = translate_sequence(&sequence, &code);
                let rscu_values = compute_rscu(&codon_counts, &code);

                // Store data for CSV output
                rscu_results.push((seq_name.clone(), rscu_values));
                codon_counts_list.push((seq_name.clone(), codon_counts));
                amino_acid_counts_list.push((seq_name.clone(), amino_acid_counts));
            }

            // Write Codon & Amino Acid Counts to CSV
            write_counts_to_csv(
                &args.output_file,
                &codon_counts_list,
                &amino_acid_counts_list,
            )
            .unwrap();

            // Write RSCU values to a single CSV file
            let rscu_filename = format!("{}_rscu.csv", args.output_file);
            write_rscu_to_csv(&rscu_filename, &rscu_results).unwrap();

            if args.compute_zscore {
                println!("Computing and saving RSCU Z-scores...");

                // Compute RSCU Z-scores
                let mean_rscu = compute_mean_rscu(&rscu_results);
                let std_rscu = compute_std_rscu(&rscu_results, &mean_rscu);
                let rscu_z_scores = compute_rscu_z_scores(&rscu_results, &mean_rscu, &std_rscu);

                // Write RSCU Z-scores to a CSV file
                let z_score_filename = format!("{}_rscu_z_scores.csv", args.output_file);
                write_z_scores_to_csv(&z_score_filename, &rscu_z_scores).unwrap();

                println!("RSCU Z-scores saved to {}", z_score_filename);
            }

            println!(
                "Results saved to:\n - {}_counts.csv\n - {}_amino_acids.csv\n - {}_rscu.csv",
                &args.output_file, &args.output_file, &args.output_file,
            );
        }
        Err(e) => eprintln!("Error reading file: {}", e),
    }
}
