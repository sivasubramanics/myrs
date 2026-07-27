use crate::error::Result;
use crate::data::fasta::{FastaReader, FastaRecord, FastaWriter};
use std::path::Path;
use crate::utils::DEFAULT_FOLD_WIDTH;

/// Generates an output filename by suffixing before the file extension.
/// e.g., "sample.fasta" -> "sample_filtered.fasta"
/// e.g., "sample.fasta.gz" -> "sample_filtered.fasta.gz"
fn generate_output_filename<P: AsRef<Path>>(input_path: P) -> String {
    let path = input_path.as_ref();
    let path_str = path.to_string_lossy();

    if path_str.ends_with(".gz") {
        let base = &path_str[..path_str.len() - 3];
        if let Some(dot_idx) = base.rfind('.') {
            format!("{}_filtered{}.gz", &base[..dot_idx], &base[dot_idx..])
        } else {
            format!("{}_filtered.gz", base)
        }
    } else if let Some(dot_idx) = path_str.rfind('.') {
        format!("{}_filtered{}", &path_str[..dot_idx], &path_str[dot_idx..])
    } else {
        format!("{}_filtered", path_str)
    }
}

/// Evaluates filtering conditions against a FASTA record.
fn passes_filters(
    record: &FastaRecord,
    min_len: Option<usize>,
    max_len: Option<usize>,
    min_gc: Option<f64>,
) -> bool {
    let len = record.len();

    if let Some(min) = min_len {
        if len < min {
            return false;
        }
    }

    if let Some(max) = max_len {
        if len > max {
            return false;
        }
    }

    if let Some(min_gc_val) = min_gc {
        if record.gc_percentage() < min_gc_val {
            return false;
        }
    }

    true
}

pub fn run(
    fname: &str,
    min_len: Option<usize>,
    max_len: Option<usize>,
    min_gc: Option<f64>,
) -> Result<()> {
    let out_fname = generate_output_filename(fname);

    // from_path and create_path return crate::error::Result
    let mut reader = FastaReader::from_path(fname)?;
    let mut writer = FastaWriter::create_path(&out_fname)?;

    let mut record = FastaRecord::default();
    let mut total_records = 0usize;
    let mut passed_records = 0usize;

    // read_next returns io::Result, but ? automatically converts io::Error -> AppError
    while reader.read_next(&mut record)? {
        total_records += 1;

        if passes_filters(&record, min_len, max_len, min_gc) {
            writer.write_record(&record, Some(DEFAULT_FOLD_WIDTH))?;
            passed_records += 1;
        }
    }

    writer.flush()?;

    println!(
        "Filtered {}/{} records written to {}",
        passed_records, total_records, out_fname
    );

    Ok(())
}