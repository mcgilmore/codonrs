//! `codonrs` is a small crate for rapidly calculating relative synonymous codon usage (RSCU)
//! values for coding DNA sequences, for analyses of codon usage bias. The crate can be used 
//! as a command-line utility with the `codonrs` command, or used in other crates via the analysis mod.
//! 
//! ## Command-line usage
//! 
//! ### Required arguments
//! 
//! **Input**: `-i`/`--input`: A multi-fasta file with the sequences to be analysed as individual
//! fasta entries. Sequences whose length is not a multiple of three will be ignored
//! 
//! **Output**: `-o`/`--output`: Prefix for the output files. Three files are output by default:
//!   - `prefix`_codon.csv: raw codon counts for each CDS.
//!   - `prefix`_amino_acids.csv: amino acid counts for each CDS, determined from the chosen translation table.
//!   - `prefix`_rscu.csv: calculated RSCU values for each CDS.
//! 
//! ### Optional arguments
//! 
//! **Translation table**: `-t`/`--table`: Integer representing NCBI translation table to be used for codon 
//! counts and RSCU calculation. See tables [here](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes).
//! Defaults to 1: the standard code.
//! 
//! 
pub mod analysis {
    use serde::{Deserialize, Serialize};
    use std::collections::{HashMap, HashSet};
    use std::error::Error;
    use std::fs::File;
    use std::io::{self, BufRead, BufReader};
    use rayon::prelude::*;
    use csv::Writer;

    static CODE_FILE: &str = include_str!("genetic_code.json");

    /// Load the provided genetic_code.json data file
    fn load_embedded_genetic_codes() -> Result<Vec<GeneticCodeJSON>, Box<dyn Error>> {
        let genetic_codes: Vec<GeneticCodeJSON> = serde_json::from_str(CODE_FILE)?;
        Ok(genetic_codes)
    }

    #[derive(Debug, Clone)]
    pub struct GeneticCode {
        id: String,
        name: String,
        codon_map: HashMap<String, String>,
    }

    /// A genetic code object, containing an id, name and translation table (codon_map)
    impl GeneticCode {
        pub fn new() -> Self {
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

    /// Get a genetic code from a GeneticCode object by a provided ID
    fn get_genetic_code_by_id<'a>(
        codes: &'a [GeneticCodeJSON],
        id: &u8,
    ) -> Option<&'a GeneticCodeJSON> {
        codes.iter().find(|code| code.id == *id)
    }

    /// Return a GeneticCode object give a translation table id (integer)
    ///
    /// # Arguments
    ///
    /// * `translation_table`: u8 representing the translation table
    ///
    /// # Returns
    ///
    /// GeneticCode object containing translation table of choice
    pub fn genetic_code_from_id(translation_table: &u8) -> GeneticCode {
        let mut code = GeneticCode::new();
        if let Ok(json_codes) = load_embedded_genetic_codes() {
            if let Some(selected_code) = get_genetic_code_by_id(&json_codes, translation_table) {
                code.id = selected_code.id.to_string();
                code.name = selected_code.name.first().cloned().unwrap_or_default();
                code.codon_map = selected_code.to_codon_map();
            } else {
                eprintln!("Error: Genetic Code ID {} not found!", translation_table);
            }
        } else {
            eprintln!("Failed to load genetic codes");
        }
        code
    }

    /// Parse a DNA sequence into codon counts.
    ///
    /// # Arguments
    ///
    /// * `sequence` - A string representing the DNA sequence.
    ///
    /// # Returns
    ///
    /// A result with a HashMap mapping codon strings to their occurrence count.
    pub fn count_codons(sequence: &str) -> HashMap<String, usize> {
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

    /// Count codons for each sequence in parallel.
    ///
    /// # Arguments
    ///
    /// * `sequences` - A vector of tuples, where each tuple contains a sequence name and the DNA sequence string.
    ///
    /// # Returns
    ///
    /// A vector of tuples, where each tuple contains the sequence name and its corresponding codon counts.
    pub fn count_codons_for_sequences(sequences: &Vec<(String, String)>) -> Vec<(String, HashMap<String, usize>)> {
        sequences.par_iter()
            .map(|(name, sequence)| (name.clone(), count_codons(sequence)))
            .collect()
    }

    /// Compute Relative Synonymous Codon Usage (RSCU) values.
    ///
    /// # Arguments
    ///
    /// * codon_counts - A hashmap with codon counts for the sequence.
    /// * code - a GeneticCode object to be used for codon translation
    ///
    /// # Returns
    ///
    /// A result with a HashMap mapping codon strings to their RSCU value
    pub fn compute_rscu(codon_counts: &HashMap<String, usize>, code: &GeneticCode) -> HashMap<String, f64> {
        use rayon::prelude::*;
        
        let codon_table = &code.codon_map;
        let mut amino_acid_totals: HashMap<&str, usize> = HashMap::new();
        let mut synonymous_codons: HashMap<&str, Vec<&str>> = HashMap::new();
        
        // Group codons by their amino acid and count occurrences
        for (codon, amino_acid) in codon_table {
            synonymous_codons
                .entry(amino_acid)
                .or_insert_with(Vec::new)
                .push(codon);
            *amino_acid_totals.entry(amino_acid).or_insert(0) += codon_counts.get(codon).copied().unwrap_or(0);
        }
        
        // Compute RSCU values in parallel
        let rscu_pairs: Vec<(String, f64)> = synonymous_codons.par_iter()
            .flat_map(|(amino_acid, codons)| {
                let total_codon_count = amino_acid_totals.get(amino_acid).copied().unwrap_or(0) as f64;
                let num_codons = codons.len() as f64;
                codons.par_iter().map(move |codon| {
                    let observed = *codon_counts.get(*codon).unwrap_or(&0) as f64;
                    let expected = total_codon_count / num_codons;
                    let rscu = if expected > 0.0 { observed / expected } else { 0.0 };
                    ((*codon).to_string(), rscu)
                })
            })
            .collect();
        
        rscu_pairs.into_iter().collect()
    }

    /// Write RSCU values to a CSV file.
    pub fn write_rscu_to_csv(
        filename: &str,
        rscu_data: &Vec<(String, HashMap<String, f64>)>,
    ) -> std::io::Result<()> {
        let mut writer = Writer::from_path(filename)?;

        // Collect all unique codons across sequences to create CSV headers
        let mut codon_set = HashSet::new();
        for (_, codon_map) in rscu_data {
            for codon in codon_map.keys() {
                codon_set.insert(codon.clone());
            }
        }

        let mut codons: Vec<String> = codon_set.into_iter().collect();
        codons.sort(); // Ensure consistent order

        // Write the header row: first column "Sequence", then one column per codon.
        let mut header = vec!["Sequence".to_string()];
        header.extend(codons.clone());
        writer.write_record(&header)?;

        // Write RSCU values for each sequence.
        for (seq_name, codon_map) in rscu_data {
            let mut record = vec![seq_name.clone()];
            for codon in &codons {
                // Format the RSCU value to ten decimal places.
                let value = codon_map.get(codon).unwrap_or(&0.0);
                record.push(format!("{:.10}", value));
            }
            writer.write_record(&record)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Compute the mean RSCU value for each codon across the genome
    pub fn compute_mean_rscu(rscu_results: &Vec<(String, HashMap<String, f64>)>) -> HashMap<String, f64> {
        let mut total_rscu: HashMap<String, f64> = HashMap::new();
        let gene_count = rscu_results.len() as f64;

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
    pub fn compute_std_rscu(
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
    pub fn compute_rscu_z_scores(
        rscu_results: &Vec<(String, HashMap<String, f64>)>,
        mean_rscu: &HashMap<String, f64>,
        std_rscu: &HashMap<String, f64>,
    ) -> Vec<(String, HashMap<String, f64>)> {
        rscu_results.par_iter()
            .map(|(gene, rscu_map)| {
                let gene_z_scores: HashMap<String, f64> = rscu_map.iter()
                    .map(|(codon, value)| {
                        let mean = mean_rscu.get(codon).unwrap_or(&0.0);
                        let std_dev = std_rscu.get(codon).unwrap_or(&1.0); // Avoid division by zero
                        let z_score = (value - mean) / std_dev;
                        (codon.clone(), z_score)
                    })
                    .collect();
                (gene.clone(), gene_z_scores)
            })
            .collect()
    }

    /// Write Z-scores to a CSV file
    pub fn write_z_scores_to_csv(
        filename: &str,
        z_scores: &Vec<(String, HashMap<String, f64>)>,
    ) -> io::Result<()> {
        let mut writer = Writer::from_path(filename)?;

        let mut codon_set = HashSet::new();
        for (_, codon_map) in z_scores {
            for codon in codon_map.keys() {
                codon_set.insert(codon.clone());
            }
        }

        let mut codons: Vec<String> = codon_set.into_iter().collect();
        codons.sort();

        // Write header.
        let mut header = vec!["Gene".to_string()];
        header.extend(codons.clone());
        writer.write_record(&header)?;

        // Write Z-scores for each gene.
        for (gene, codon_map) in z_scores {
            let mut record = vec![gene.clone()];
            for codon in &codons {
                let z_score = codon_map.get(codon).unwrap_or(&0.0);
                record.push(format!("{:.10}", z_score));
            }
            writer.write_record(&record)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Count the amino acids for a CDS
    ///
    /// # Arguments
    ///
    /// * sequence: a str of the sequence to be analysed
    /// * code: a GeneticCode object to be used for codon translation
    ///
    /// # Returns
    ///
    /// A result with a HashMap mapping amino acids to their count
    pub fn count_amino_acids(sequence: &str, code: &GeneticCode) -> HashMap<String, usize> {
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

    /// Write codon counts to a CSV file
    ///
    /// # Arguments
    ///
    /// * filename_prefix: str to be used as prefix for output files
    /// * codon_data: The codon counts to be written
    /// 
    /// # Returns
    ///
    /// A result with an output file: `prefix`_codon.csv
    pub fn write_codon_counts_to_csv(
        filename_prefix: &str,
        codon_data: &Vec<(String, HashMap<String, usize>)>,
    ) -> std::io::Result<()> {

        let codon_filename = format!("{}_codon.csv", filename_prefix);
        let mut codon_writer = Writer::from_path(codon_filename)?;
    
        // Collect all unique codons across sequences.
        let mut codon_set = HashSet::new();
        for (_, codon_map) in codon_data {
            for codon in codon_map.keys() {
                codon_set.insert(codon.clone());
            }
        }
        let mut codons: Vec<String> = codon_set.into_iter().collect();
        codons.sort(); // Ensure consistent order
    
        // Write the header row: "Sequence", then one column per codon.
        let mut header = vec!["Sequence".to_string()];
        header.extend(codons.clone());
        codon_writer.write_record(&header)?;
    
        // Write codon counts for each sequence.
        for (seq_name, codon_map) in codon_data {
            let mut record = vec![seq_name.clone()];
            for codon in &codons {
                let count = codon_map.get(codon).unwrap_or(&0);
                record.push(count.to_string());
            }
            codon_writer.write_record(&record)?;
        }
        codon_writer.flush()?;
        Ok(())
    }

    /// Write amino acid counts to a CSV file
    ///
    /// # Arguments
    ///
    /// * filename_prefix: str to be used as prefix for output files
    /// * amino_acid_data: The amino acid counts to be written
    /// 
    /// # Returns
    ///
    /// A result with an output files: `prefix`_amino_acids.csv for amino acids
    pub fn write_amino_acid_counts_to_csv(
        filename_prefix: &str,
        amino_acid_data: &Vec<(String, HashMap<String, usize>)>,
    ) -> std::io::Result<()> {
        let amino_filename = format!("{}_amino_acids.csv", filename_prefix);
        let mut amino_writer = Writer::from_path(amino_filename)?;
    
        // Collect all unique amino acids across sequences.
        let mut amino_acid_set = HashSet::new();
        for (_, amino_map) in amino_acid_data {
            for amino in amino_map.keys() {
                amino_acid_set.insert(amino.clone());
            }
        }
        let mut amino_acids: Vec<String> = amino_acid_set.into_iter().collect();
        amino_acids.sort(); // Ensure consistent order
    
        // Write the header row: "Sequence", then one column per amino acid.
        let mut header = vec!["Sequence".to_string()];
        header.extend(amino_acids.clone());
        amino_writer.write_record(&header)?;
    
        // Write amino acid counts for each sequence.
        for (seq_name, amino_map) in amino_acid_data {
            let mut record = vec![seq_name.clone()];
            for amino in &amino_acids {
                let count = amino_map.get(amino).unwrap_or(&0);
                record.push(count.to_string());
            }
            amino_writer.write_record(&record)?;
        }
        amino_writer.flush()?;
        Ok(())
    }


    /// Read sequences from a multi-FASTA file
    ///
    /// # Arguments
    ///
    /// * filename_prefix: str to be used as prefix for output files
    /// * codon_data: The codon counts to be written
    /// * amino_acid_data: The amino acid counts to be written
    /// 
    /// # Returns
    ///
    /// A result with two output files: `prefix`_codon.csv for codons and `prefix`_amino_acids.csv for amino acids
    pub fn read_sequences_from_fasta(filename: &str) -> io::Result<Vec<(String, String)>> {
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
}
