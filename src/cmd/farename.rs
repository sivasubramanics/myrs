use crate::data::fasta::{FastaReader, FastaRecord, FastaWriter};
use crate::error::{AppError, Result};
use crate::utils::DEFAULT_FOLD_WIDTH;
use log::{debug, info, warn};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

fn generate_output_filename<P: AsRef<Path>>(input_path: P) -> String {
    let path = input_path.as_ref();
    let path_str = path.to_string_lossy();

    if path_str.ends_with(".gz") {
        let base = &path_str[..path_str.len() - 3];
        if let Some(dot_idx) = base.rfind('.') {
            format!("{}.renamed{}.gz", &base[..dot_idx], &base[dot_idx..])
        } else {
            format!("{}.renamed.gz", base)
        }
    } else if let Some(dot_idx) = path_str.rfind('.') {
        format!("{}.renamed{}", &path_str[..dot_idx], &path_str[dot_idx..])
    } else {
        format!("{}.renamed", path_str)
    }
}

/// Helper function to parse a 2-column TSV mapping file into a HashMap (old_id -> new_id)
fn load_name_map<P: AsRef<Path>>(map_path: P) -> Result<HashMap<String, String>> {
    let file = File::open(&map_path)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 2 {
            warn!(
                "line {} in mapping file '{:?}' has fewer than 2 tab-separated columns; skipping",
                line_idx + 1,
                map_path.as_ref()
            );
            continue;
        }

        let old_name = fields[0].trim().to_string();
        let new_name = fields[1].trim().to_string();

        if !old_name.is_empty() && !new_name.is_empty() {
            map.insert(old_name, new_name);
        }
    }

    Ok(map)
}

pub fn run(
    fname: &str,
    map_file: &str,
    output: Option<String>,
    keep_old_id: bool,
    fold_width: Option<usize>,
) -> Result<()> {
    let name_map = load_name_map(map_file)?;
    if name_map.is_empty() {
        warn!("mapping file {} is empty or contains no valid pairs; records will remain unchanged", map_file);
    }

    let out_fname = output.unwrap_or_else(|| generate_output_filename(fname));
    let wrap_width = fold_width.or(Some(DEFAULT_FOLD_WIDTH));

    info!("renaming sequence headers in input file: {}", fname);
    debug!(
        "rename config - map_file: {}, keep_old_id: {}, output: {}",
        map_file, keep_old_id, out_fname
    );

    let mut reader = FastaReader::from_path(fname)?;
    let mut writer = FastaWriter::create_path(&out_fname)?;

    let mut record = FastaRecord::default();
    let mut total_records = 0usize;
    let mut renamed_records = 0usize;

    while reader.read_next(&mut record)? {
        total_records += 1;

        if let Some(new_name) = name_map.get(&record.name) {
            let old_name = std::mem::replace(&mut record.name, new_name.clone());

            if keep_old_id {
                let tag = format!("old_id={}", old_name);
                if record.description.trim().is_empty() {
                    record.description = tag;
                } else {
                    record.description.push(' ');
                    record.description.push_str(&tag);
                }
            }
            renamed_records += 1;
        }

        writer.write_record(&record, wrap_width)?;
    }

    writer.flush()?;

    if total_records == 0 {
        warn!("input file {} contains zero records.", fname);
    } else {
        info!(
            "processed {} sequence(s). Renamed {}/{} records written to {}",
            total_records, renamed_records, total_records, out_fname
        );
    }

    Ok(())
}