use crate::data::fasta::FastaReader;
use crate::error::{AppError, Result};
use crate::utils::defaults::DEFAULT_COLUMN_WIDTH;
use crate::utils::helpers::{create_file, num_to_str, nx};
use log::{debug, info, warn};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn run(fname: &str) -> Result<()> {
    info!("Starting summary analysis on file: {}", fname);

    let reader = FastaReader::from_path(fname)?;

    let mut total_count: usize = 0;
    let mut total_length: usize = 0;
    let mut max_length: usize = 0;
    let mut min_length: usize = usize::MAX;

    let mut non_atgc: usize = 0;
    let mut lengths: Vec<usize> = Vec::new();

    let mut ge100 = 0;
    let mut ge200 = 0;
    let mut ge500 = 0;
    let mut ge1k = 0;
    let mut ge2k = 0;
    let mut ge10k = 0;
    let mut ge1m = 0;

    for record in reader {
        let record = record?;
        let len = record.len();

        total_count += 1;
        total_length += len;
        max_length = max_length.max(len);
        min_length = min_length.min(len);

        lengths.push(len);

        match len {
            l if l >= 1_000_000 => ge1m += 1,
            l if l >= 10_000 => ge10k += 1,
            l if l >= 2_000 => ge2k += 1,
            l if l >= 1_000 => ge1k += 1,
            l if l >= 500 => ge500 += 1,
            l if l >= 200 => ge200 += 1,
            l if l >= 100 => ge100 += 1,
            _ => {}
        }

        non_atgc += record.num_n();
    }

    if total_count == 0 {
        warn!("Input file '{}' contained 0 records.", fname);
        min_length = 0;
    } else {
        debug!("Parsed {} records. Sorting contig lengths...", total_count);
    }

    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let n50 = nx(&lengths, total_length, 0.50);
    let n90 = nx(&lengths, total_length, 0.90);

    let avg = if total_count > 0 {
        total_length / total_count
    } else {
        0
    };

    let pct = if total_length > 0 {
        (non_atgc as f64 * 100.0) / total_length as f64
    } else {
        0.0
    };

    // --- 1. Formatted Summary Block via Logger ---
    // --- 1. Bioinformatics-Standard Summary Block ---
    info!("====================== SUMMARY ======================");
    info!("{:<width$} : {}", "Input file", fname, width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Total sequences", num_to_str(total_count as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {} bp", "Total length", num_to_str(total_length as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {} bp", "Min length", num_to_str(min_length as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {} bp", "Avg length", num_to_str(avg as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {} bp", "Max length", num_to_str(max_length as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("------------------ COMPOSITION ------------------");
    info!("{:<width$} : {}", "Non-ATGC bases (N)", num_to_str(non_atgc as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {:.2}%", "Non-ATGC (%)", pct, width = DEFAULT_COLUMN_WIDTH);
    info!("------------------- CONTIGUITY ------------------");
    info!("{:<width$} : {} bp", "N50", num_to_str(n50 as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {} bp", "N90", num_to_str(n90 as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("------------------- LENGTH BINS -----------------");
    info!("{:<width$} : {}", "Count >= 100 bp", num_to_str(ge100 as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Count >= 200 bp", num_to_str(ge200 as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Count >= 500 bp", num_to_str(ge500 as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Count >= 1 Kb", num_to_str(ge1k as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Count >= 2 Kb", num_to_str(ge2k as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Count >= 10 Kb", num_to_str(ge10k as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("{:<width$} : {}", "Count >= 1 Mb", num_to_str(ge1m as u64), width = DEFAULT_COLUMN_WIDTH);
    info!("=================================================");

    // --- 2. Write TSV Output File (`input.fasta.summary.tsv`) ---
    let tsv_path = format!("{}.summary.tsv", fname);
    let file = create_file(&tsv_path)?;

    let mut writer = BufWriter::new(file);

    writeln!(writer, "input_file\t{}", fname)?;
    writeln!(writer, "num_sequences\t{}", total_count)?;
    writeln!(writer, "total_length\t{}", total_length)?;
    writeln!(writer, "min_length\t{}", min_length)?;
    writeln!(writer, "avg_length\t{}", avg)?;
    writeln!(writer, "max_length\t{}", max_length)?;
    writeln!(writer, "non_atgc_bases\t{}", non_atgc)?;
    writeln!(writer, "non_atgc_pct\t{:.2}", pct)?;
    writeln!(writer, "n50\t{}", n50)?;
    writeln!(writer, "n90\t{}", n90)?;
    writeln!(writer, "count_ge_100bp\t{}", ge100)?;
    writeln!(writer, "count_ge_200bp\t{}", ge200)?;
    writeln!(writer, "count_ge_500bp\t{}", ge500)?;
    writeln!(writer, "count_ge_1kb\t{}", ge1k)?;
    writeln!(writer, "count_ge_2kb\t{}", ge2k)?;
    writeln!(writer, "count_ge_10kb\t{}", ge10k)?;
    writeln!(writer, "count_ge_1mb\t{}", ge1m)?;

    writer.flush()?;

    info!("Saved summary metrics to TSV report: {}", tsv_path);

    Ok(())
}