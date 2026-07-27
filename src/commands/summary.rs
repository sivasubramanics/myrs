use crate::data::fasta::FastaReader;
use crate::utils::helpers::{num_to_str, nx};
use crate::utils::defaults::DEFAULT_COLUMN_WIDTH;


pub fn run(fname: &str) -> Result<(), Box<dyn std::error::Error>> {
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

    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let n50 = nx(&lengths, total_length, 0.50);
    let n90 = nx(&lengths, total_length, 0.90);

    // Adjust COL_WIDTH to whatever width best fits your terminal/needs

    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "Total Contigs:", num_to_str(total_count as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "Maximum Contig Length:", num_to_str(max_length as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "Minimum Contig Length:", num_to_str(min_length as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "Total Contig Length:", num_to_str(total_length as u64));

    if total_count > 0 {
        println!(
            "{:<DEFAULT_COLUMN_WIDTH$}{:.2}",
            "Average Contig Length:",
            num_to_str(total_length as u64 / total_count as u64)
        );
    }

    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "Total Number of nonATGC:", num_to_str(non_atgc as u64));

    let pct = if total_length > 0 {
        (non_atgc as f64 * 100.0) / total_length as f64
    } else {
        0.0
    };

    println!("{:<DEFAULT_COLUMN_WIDTH$}{:.2}%", "Percentage nonATGC:", pct);

    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=100bp:", num_to_str(ge100 as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=200bp:", num_to_str(ge200 as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=500bp:", num_to_str(ge500 as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=1kb:", num_to_str(ge1k as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=2kb:", num_to_str(ge2k as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=10kb:", num_to_str(ge10k as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "No. Contigs >=1mb:", num_to_str(ge1m as u64));

    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "N50 Value:", num_to_str(n50 as u64));
    println!("{:<DEFAULT_COLUMN_WIDTH$}{}", "N90 Value:", num_to_str(n90 as u64));

    Ok(())
}

