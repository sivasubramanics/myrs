use crate::data::fasta::{FastaRecord, FastaWriter};
use crate::error::{AppError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct FaidxRecord {
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
    pub description: String,
}

pub struct FaidxIndex {
    records: HashMap<String, FaidxRecord>,
    fasta_path: PathBuf,
}

impl FaidxIndex {
    /// Load `.faidx` index file into memory and associate it with its source FASTA file.
    pub fn from_path<P1: AsRef<Path>, P2: AsRef<Path>>(fai_path: P1, fasta_path: P2) -> Result<Self> {
        let file = File::open(fai_path.as_ref())?;
        let reader = BufReader::new(file);
        let mut records = HashMap::new();

        for line in reader.lines() {
            let l = line?;
            if l.trim().is_empty() {
                continue;
            }

            let parts: Vec<&str> = l.split('\t').collect();
            if parts.len() < 6 {
                return Err(AppError::InvalidInput("malformed .faidx entry".into()));
            }

            let description = if parts.len() > 5 && parts[5] != "." {
                parts[5].to_string()
            } else {
                String::new()
            };

            let record = FaidxRecord {
                name: parts[0].to_string(),
                length: parts[1].parse().map_err(|_| AppError::InvalidInput("invalid length".into()))?,
                offset: parts[2].parse().map_err(|_| AppError::InvalidInput("invalid offset".into()))?,
                line_bases: parts[3].parse().map_err(|_| AppError::InvalidInput("invalid line_bases".into()))?,
                line_width: parts[4].parse().map_err(|_| AppError::InvalidInput("invalid line_width".into()))?,
                description,
            };

            records.insert(record.name.clone(), record);
        }

        Ok(Self {
            records,
            fasta_path: fasta_path.as_ref().to_path_buf(),
        })
    }

    /// Retrieve the full sequence for a given sequence ID as a `FastaRecord`.
    pub fn fetch(&self, id: &str) -> Result<FastaRecord> {
        let meta = self.records.get(id).ok_or_else(|| {
            AppError::InvalidInput(format!("record '{}' not found in index", id))
        })?;

        self.read_range(meta, 0, meta.length as usize)
    }

    /// Retrieve a region slice (`start` index with a given `length`) as a `FastaRecord`.
    /// `start` is 0-based. If `start + length` exceeds total sequence length, it gets safely capped.
    pub fn fetch_range(&self, id: &str, start: usize, length: usize) -> Result<FastaRecord> {
        let meta = self.records.get(id).ok_or_else(|| {
            AppError::InvalidInput(format!("record '{}' not found in index", id))
        })?;

        if start >= meta.length as usize {
            return Err(AppError::InvalidInput(format!(
                "start coordinate {} out of bounds for sequence '{}' (len: {})",
                start, id, meta.length
            )));
        }

        // Cap length if it extends past sequence boundary
        let actual_len = length.min(meta.length as usize - start);

        self.read_range(meta, start, actual_len)
    }

    /// Core helper to seek and slice sequence data directly from disk.
    fn read_range(&self, meta: &FaidxRecord, start: usize, length: usize) -> Result<FastaRecord> {
        let start = start as u64;
        let length = length as u64;

        if length == 0 {
            return Ok(FastaRecord {
                id: 0,
                name: format!("{}:{}-{}", meta.name, start, start),
                description: String::new(),
                sequence: Vec::new(),
            });
        }

        let newline_bytes = meta.line_width - meta.line_bases;

        // Calculate exact byte offset on disk corresponding to `start`
        let lines_before = start / meta.line_bases;
        let bases_on_last_line = start % meta.line_bases;
        let byte_offset = meta.offset + (lines_before * meta.line_width) + bases_on_last_line;

        // Calculate total raw bytes to read (bases + intervening line endings)
        let end_base = start + length;
        let start_line = start / meta.line_bases;
        let end_line = (end_base - 1) / meta.line_bases;
        let num_line_switches = end_line - start_line;
        let raw_bytes_to_read = length + (num_line_switches * newline_bytes);

        // Perform single disk seek & read
        let mut file = File::open(&self.fasta_path)?;
        file.seek(SeekFrom::Start(byte_offset))?;

        let mut raw_buf = vec![0u8; raw_bytes_to_read as usize];
        file.read_exact(&mut raw_buf)?;

        // Filter out line ending bytes (`\n` and `\r`)
        let sequence: Vec<u8> = raw_buf
            .into_iter()
            .filter(|&b| b != b'\n' && b != b'\r')
            .collect();

        // Format name tag like standard samtools output (e.g. chr01:100-200) if a subset range was requested
        let record_name = if length == meta.length && start == 0 {
            meta.name.clone()
        } else {
            format!("{}:{}-{}", meta.name, start + 1, start + length)
        };

        Ok(FastaRecord {
            id: 0,
            name: record_name,
            description: meta.description.clone(),
            sequence,
        })
    }
}