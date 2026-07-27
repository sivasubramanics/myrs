use crate::error::{AppError, Result};
use log::{error, info, warn};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use crate::utils::helpers::{create_file, open_file, is_compressed};

pub fn run(fname: &str) -> Result<()> {
    info!("creating fai file '{}'", fname);
    // if the input file compressed, throw error saying we can't index compressed files

    if is_compressed(fname) {
        error!("Input file '{}' is compressed. Please provide an uncompressed FASTA file for indexing.", fname);
        return Err(AppError::InvalidInput("Compressed file not supported".into()));
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
            // End of file — flush the last active sequence
            if let Some(name) = seq_name.take() {
                let desc = seq_desc.take().unwrap_or_default();
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_width, desc
                )?;
            }
            break;
        }

        if line_buf.starts_with(b">") {
            // If we encounter a new header, write the prior record's index entry
            if let Some(name) = seq_name.take() {
                let desc = seq_desc.take().unwrap_or_default();
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_width, desc
                )?;
            }

            // Parse header string
            let header_str = String::from_utf8_lossy(&line_buf[1..]).trim_end().to_string();
            let mut parts = header_str.splitn(2, |c: char| c.is_whitespace());

            seq_name = Some(parts.next().unwrap_or("").to_string());
            seq_desc = Some(parts.next().unwrap_or("").to_string());

            seq_len = 0;
            seq_offset = current_offset + bytes_read as u64;
            line_bases = 0;
            line_width = 0;
            is_first_seq_line = true;
        } else if seq_name.is_some() {
            // Sequence lines logic
            let raw_len = line_buf.len() as u64;
            if raw_len > 0 {
                // Determine newline length (\n vs \r\n)
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