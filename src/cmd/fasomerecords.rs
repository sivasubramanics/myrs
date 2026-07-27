use crate::cmd::fafai;
use crate::data::faidx::FaidxIndex;
use crate::data::fasta::{FastaReader, FastaRecord, FastaWriter};
use crate::error::{AppError, Result};
use log::{info, warn};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Threshold where sequential streaming becomes faster/more efficient than thousands of random seeks
const SEQUENTIAL_THRESHOLD: usize = 1000;

pub fn run(
    fname: &str,
    names_file: Option<&str>,
    names_list: Option<Vec<String>>,
    output: Option<String>,
    fold_width: Option<usize>,
) -> Result<()> {
    let targets = load_target_names(names_file, names_list)?;
    if targets.is_empty() {
        return Err(AppError::InvalidInput("no sequence names provided to extract".into()));
    }

    let fasta_path = Path::new(fname);
    let fai_path = PathBuf::from(format!("{}.faidx", fname));
    let out_path = match output {
        Some(path_str) => PathBuf::from(path_str),
        None => derive_output_path(fasta_path)?,
    };

    let mut writer = FastaWriter::create_path(&out_path)?;

    if targets.len() > SEQUENTIAL_THRESHOLD {
        info!(
            "target count ({}) exceeds threshold ({}); using sequential streaming via FastaReader...",
            targets.len(),
            SEQUENTIAL_THRESHOLD
        );
        extract_sequential(fasta_path, &targets, &mut writer, fold_width)?;
    } else {
        if !fai_path.exists() {
            info!("index file {:?} not found. Generating .faidx index...", fai_path);
            fafai::run(fname)?;
        }

        info!("extracting {} records via faidx random access...", targets.len());
        extract_indexed(&fai_path, fasta_path, &targets, &mut writer, fold_width)?;
    }

    writer.flush()?;
    info!("successfully wrote extracted records to {:?}", out_path);
    Ok(())
}

/// Fast O(1) indexed extraction using faidx seeks
fn extract_indexed(
    fai_path: &Path,
    fasta_path: &Path,
    targets: &[String],
    writer: &mut FastaWriter<Box<dyn std::io::Write>>,
    fold_width: Option<usize>,
) -> Result<()> {
    let index = FaidxIndex::from_path(fai_path, fasta_path)?;
    let mut found_count = 0;

    for name in targets {
        match index.fetch(name) {
            Ok(record) => {
                writer.write_record(&record, fold_width)?;
                found_count += 1;
            }
            Err(_) => {
                warn!("sequence record '{}' not found in faidx index", name);
            }
        }
    }

    info!("extracted {}/{} requested records", found_count, targets.len());
    Ok(())
}

/// Single-pass sequential streaming reader for large query sets
fn extract_sequential(
    fasta_path: &Path,
    targets: &[String],
    writer: &mut FastaWriter<Box<dyn std::io::Write>>,
    fold_width: Option<usize>,
) -> Result<()> {
    let mut reader = FastaReader::from_path(fasta_path)?;
    let mut record = FastaRecord::default();

    let mut remaining_targets: HashSet<&str> = targets.iter().map(|s| s.as_str()).collect();

    let mut matched_records: std::collections::HashMap<String, FastaRecord> =
        std::collections::HashMap::new();

    while reader.read_next(&mut record)? {
        if remaining_targets.contains(record.name.as_str()) {
            remaining_targets.remove(record.name.as_str());
            matched_records.insert(record.name.clone(), record.clone());

            if remaining_targets.is_empty() {
                break;
            }
        }
    }

    let mut found_count = 0;
    for name in targets {
        if let Some(record) = matched_records.get(name) {
            writer.write_record(record, fold_width)?;
            found_count += 1;
        } else {
            warn!("sequence record '{}' not found in FASTA", name);
        }
    }

    info!("extracted {}/{} requested records", found_count, targets.len());
    Ok(())
}

/// Helper function to parse target names from file or command-line list
/// Reads sequence names while preserving input order and skipping duplicates
fn load_target_names(
    file_path: Option<&str>,
    inline_list: Option<Vec<String>>,
) -> Result<Vec<String>> {
    let mut ordered_names = Vec::new();
    let mut seen = HashSet::new();

    let mut add_name = |name: String| {
        if seen.insert(name.clone()) {
            ordered_names.push(name);
        }
    };

    if let Some(list) = inline_list {
        for name in list {
            let trimmed = name.trim().to_string();
            if !trimmed.is_empty() {
                add_name(trimmed);
            }
        }
    }

    if let Some(path_str) = file_path {
        let file = File::open(path_str)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                add_name(trimmed.to_string());
            }
        }
    }

    Ok(ordered_names)
}

fn derive_output_path(input_path: &Path) -> Result<PathBuf> {
    let parent = input_path.parent().unwrap_or_else(|| Path::new(""));
    let file_name = input_path
        .file_name()
        .and_then(|f| f.to_str())
        .ok_or_else(|| AppError::InvalidInput(format!("invalid input path '{:?}'", input_path)))?;

    let new_filename = if let Some((stem, ext)) = file_name.rsplit_once('.') {
        format!("{}.some.{}", stem, ext)
    } else {
        format!("{}.some", file_name)
    };

    Ok(parent.join(new_filename))
}