use codonrs::analysis::*;
use std::collections::HashMap;
use rayon::prelude::*;

#[test]
fn sample_file() {
    let filename: &str = "tests/cds_adapted.fasta";
    let sequences = read_sequences_from_fasta(&filename).unwrap();
    let translation_table: u8 = 2;
    let code = genetic_code_from_id(&translation_table);
    
    let results: Vec<(String, HashMap<String, usize>, HashMap<String, usize>, HashMap<String, f64>)> =
        count_codons_for_sequences(&sequences)
            .into_par_iter()
            .filter_map(|(seq_name, codon_counts)| {
                if codon_counts.is_empty() {
                    eprintln!(
                        "Warning: Skipping sequence '{}' (not a multiple of 3)",
                        seq_name
                    );
                    None
                } else {
                    // Retrieve the corresponding sequence from the original vector
                    let sequence = sequences
                        .iter()
                        .find(|(name, _)| name.as_str() == seq_name.as_str())
                        .map(|(_, seq)| seq)
                        .expect("Sequence not found");
                    
                    // Compute additional metrics for this sequence.
                    let amino_acid_counts = count_amino_acids(sequence, &code);
                    let rscu_values = compute_rscu(&codon_counts, &code);
                    
                    Some((seq_name, codon_counts, amino_acid_counts, rscu_values))
                }
            })
            .collect();

    // Explicitly annotate the types for output vectors.
    let mut codon_counts_list: Vec<(String, HashMap<String, usize>)> = Vec::new();
    let mut amino_acid_counts_list: Vec<(String, HashMap<String, usize>)> = Vec::new();
    let mut rscu_results: Vec<(String, HashMap<String, f64>)> = Vec::new();

    // Split the results into separate vectors for further analysis or CSV output.
    for (seq_name, codon_counts, amino_acid_counts, rscu_values) in results.into_iter() {
        codon_counts_list.push((seq_name.clone(), codon_counts));
        amino_acid_counts_list.push((seq_name.clone(), amino_acid_counts));
        rscu_results.push((seq_name, rscu_values));
    }
    
    // Now, analyze the output and only pass the test if it matches the expected results.
    // For example, let's assert that the number of genes is as expected.
    let expected_number_of_genes = 5355; // change to your expected number
    assert_eq!(rscu_results.len(), expected_number_of_genes, "Number of genes mismatch");

    // Define a tolerance for floating-point comparisons.
    let tolerance = 1e-6;
    
    // Example: Check that for gene "gene1", the RSCU value for codon "ATG" is approximately 1.0.
    for (gene, rscu_map) in &rscu_results {
        if gene == "Atu0025" {
            // Ensure the codon exists in the map.
            let atc_rscu = rscu_map.get("ATC").expect("Codon ATC not found in gene1");
            let expected_atc_rscu = 2.294118; // change this to your expected value
            assert!(
                (atc_rscu - expected_atc_rscu).abs() < tolerance,
                "RSCU for ATC in {} is not within tolerance: expected {}, got {}",
                gene,
                expected_atc_rscu,
                atc_rscu
            );
        }
    }
    
    // Additional assertions can be added as needed.
}