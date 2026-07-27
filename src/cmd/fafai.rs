use crate::error::{AppError, Result};
use log::{error, info, warn};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use crate::data::fasta::FastaReader;
use crate::utils::helpers::{create_file, open_file, is_compressed, has_valid_extension};
use crate::utils::VALID_FASTA_EXTENSIONS;

pub fn run(fname: &str) -> Result<()> {
    info!("creating faidx file '{}'", fname);

    if !has_valid_extension(fname, VALID_FASTA_EXTENSIONS, false) {
        error!("input file '{}' does not have a valid FASTA extension. Supported extensions: .fasta, .fa, .fna, .faa, .ffn, .frn (not compressed)", fname);
        return Err(AppError::InvalidInput("invalid file extension".into()));
    }

    let file = open_file(fname)?;
    let mut reader = BufReader::new(file);

    let out_path = format!("{}.faidx", fname);
    let out_file = create_file(&out_path)?;
    let mut writer = BufWriter::new(out_file);

    let mut current_offset: u64 = 0;
    let mut line_buf = Vec::new();

    let mut seq_name: Option<String> = None;
    let mut seq_desc: Option<String> = None;
    let mut seq_len: u64 = 0;
    let mut seq_offset: u64 = 0;
    let mut line_bases: u64 = 0;
    let mut line_width: u64 = 0;
    let mut is_first_seq_line = true;

    loop {
        line_buf.clear();
        let bytes_read = reader.read_until(b'\n', &mut line_buf).map_err(|e| AppError::FileIo {
            path: Path::new(fname).to_path_buf(),
            source: e,
        })?;

        if bytes_read == 0 {
            // End of file — flush the last sequence
            if let Some(name) = seq_name.take() {
                let desc = seq_desc
                    .take()
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| ".".to_string());

                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_width, desc
                )?;
            }
            break;
        }

        if line_buf.starts_with(b">") {
            // New header found — flush prior record
            if let Some(name) = seq_name.take() {
                let desc = seq_desc
                    .take()
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| ".".to_string());

                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_width, desc
                )?;
            }

            let header_str = String::from_utf8_lossy(&line_buf[1..]).trim_end().to_string();
            let mut parts = header_str.splitn(2, |c: char| c.is_whitespace());

            seq_name = parts.next().map(|s| s.to_string());
            seq_desc = parts.next().map(|s| s.trim_start().replace('\t', " "));

            seq_len = 0;
            seq_offset = current_offset + bytes_read as u64;
            line_bases = 0;
            line_width = 0;
            is_first_seq_line = true;
        } else if seq_name.is_some() {
            let raw_len = line_buf.len() as u64;
            if raw_len > 0 {
                let newline_len = if line_buf.ends_with(b"\r\n") {
                    2
                } else if line_buf.ends_with(b"\n") {
                    1
                } else {
                    0
                };

                let bases_in_line = raw_len - newline_len;

                if bases_in_line > 0 {
                    seq_len += bases_in_line;

                    if is_first_seq_line {
                        line_bases = bases_in_line;
                        line_width = raw_len;
                        is_first_seq_line = false;
                    }
                }
            }
        }

        current_offset += bytes_read as u64;
    }

    writer.flush()?;
    info!("faidx file successfully written to {}", out_path);

    Ok(())
}