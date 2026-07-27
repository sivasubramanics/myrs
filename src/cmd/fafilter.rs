use crate::data::fasta::{FastaReader, FastaRecord, FastaWriter};
use crate::error::Result;
use crate::utils::DEFAULT_FOLD_WIDTH;
use log::{debug, info, warn};
use std::path::Path;

fn generate_output_filename<P: AsRef<Path>>(input_path: P) -> String {
    let path = input_path.as_ref();
    let path_str = path.to_string_lossy();

    if path_str.ends_with(".gz") {
        let base = &path_str[..path_str.len() - 3];
        if let Some(dot_idx) = base.rfind('.') {
            format!("{}.filtered{}.gz", &base[..dot_idx], &base[dot_idx..])
        } else {
            format!("{}.filtered.gz", base)
        }
    } else if let Some(dot_idx) = path_str.rfind('.') {
        format!("{}.filtered{}", &path_str[..dot_idx], &path_str[dot_idx..])
    } else {
        format!("{}.filtered", path_str)
    }
}

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

    info!("filtering on input file: {}", fname);
    debug!(
        "filter criteria applied - min_len: {:?}, max_len: {:?}, min_gc: {:?}",
        min_len, max_len, min_gc
    );

    let mut reader = FastaReader::from_path(fname)?;
    let mut writer = FastaWriter::create_path(&out_fname)?;

    let mut record = FastaRecord::default();
    let mut total_records = 0usize;
    let mut passed_records = 0usize;

    while reader.read_next(&mut record)? {
        total_records += 1;

        if passes_filters(&record, min_len, max_len, min_gc) {
            writer.write_record(&record, Some(DEFAULT_FOLD_WIDTH))?;
            passed_records += 1;
        }
    }

    writer.flush()?;

    let removed_records = total_records - passed_records;

    if total_records == 0 {
        warn!("input file {} contains zero records.", fname);
    } else {
        info!(
            "removed {} sequence(s) from file. Kept {}/{} records written to {}",
            removed_records, passed_records, total_records, out_fname
        );
    }

    Ok(())
}