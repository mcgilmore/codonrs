use clap::Parser;
use tqdm::Iter;
use codonrs::analysis;

/// Command-line arguments using Clap
#[derive(Parser)]
#[command(name = "codonrs")]
#[command(version = "0.2.1")]
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

fn main() {
    let args = Cli::parse();

    let mut rscu_results = Vec::new();
    let mut codon_counts_list = Vec::new();
    let mut amino_acid_counts_list = Vec::new();

    let code = analysis::genetic_code_from_id(&args.translation_table);
    
    match analysis::read_sequences_from_fasta(&args.input_file) {
        Ok(sequences) => {
            print!("Calculating codon counts and RSCU values for {} sequences...\n", sequences.len());
            for (seq_name, sequence) in sequences.iter().tqdm() {
                let codon_counts = analysis::count_codons(&sequence);

                if codon_counts.is_empty() {
                    eprintln!(
                        "Warning: Skipping sequence '{}' (not a multiple of 3)",
                        seq_name
                    );
                    continue;
                }

                let amino_acid_counts = analysis::count_amino_acids(&sequence, &code);
                let rscu_values = analysis::compute_rscu(&codon_counts, &code);

                // Store data for CSV output
                rscu_results.push((seq_name.clone(), rscu_values));
                codon_counts_list.push((seq_name.clone(), codon_counts));
                amino_acid_counts_list.push((seq_name.clone(), amino_acid_counts));
            }

            print!("Writing results to CSV...\n");
            // Write Codon & Amino Acid Counts to CSV
            analysis::write_counts_to_csv(
                &args.output_file,
                &codon_counts_list,
                &amino_acid_counts_list,
            )
            .unwrap();

            // Write RSCU values to a single CSV file
            let rscu_filename = format!("{}_rscu.csv", args.output_file);
            analysis::write_rscu_to_csv(&rscu_filename, &rscu_results).unwrap();

            if args.compute_zscore {
                println!("Computing and saving RSCU Z-scores...");

                // Compute RSCU Z-scores
                let mean_rscu = analysis::compute_mean_rscu(&rscu_results);
                let std_rscu = analysis::compute_std_rscu(&rscu_results, &mean_rscu);
                let rscu_z_scores = analysis::compute_rscu_z_scores(&rscu_results, &mean_rscu, &std_rscu);

                // Write RSCU Z-scores to a CSV file
                let z_score_filename = format!("{}_rscu_z_scores.csv", args.output_file);
                analysis::write_z_scores_to_csv(&z_score_filename, &rscu_z_scores).unwrap();

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
