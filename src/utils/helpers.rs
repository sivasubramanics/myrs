use std::fs::File;
use std::path::Path;
use crate::data::fasta::FastaRecord;
use crate::error::{AppError, Result};

/// Trims leading and trailing ASCII whitespace from a byte slice in-place without allocation.
#[inline]
pub fn trim_ascii_whitespace(mut bytes: &[u8]) -> &[u8] {
    while let Some((first, rest)) = bytes.split_first() {
        if first.is_ascii_whitespace() {
            bytes = rest;
        } else {
            break;
        }
    }
    while let Some((last, rest)) = bytes.split_last() {
        if last.is_ascii_whitespace() {
            bytes = rest;
        } else {
            break;
        }
    }
    bytes
}

/// Generates a modified output path by appending a suffix before extension(s).
pub fn append_suffix_to_path<P: AsRef<Path>>(input_path: P, suffix: &str) -> String {
    let path = input_path.as_ref();
    let path_str = path.to_string_lossy();

    if path_str.ends_with(".gz") {
        let base = &path_str[..path_str.len() - 3];
        if let Some(dot_idx) = base.rfind('.') {
            format!("{}{}{}.gz", &base[..dot_idx], suffix, &base[dot_idx..])
        } else {
            format!("{}{}.gz", base, suffix)
        }
    } else if let Some(dot_idx) = path_str.rfind('.') {
        format!("{}{}{}", &path_str[..dot_idx], suffix, &path_str[dot_idx..])
    } else {
        format!("{}{}", path_str, suffix)
    }
}

/// Generic checker for file extensions with optional trailing `.gz`.
///
/// Accepts a path and a list of valid extensions (e.g. `["fasta", "fa", "fna"]`).
/// Handles lowercase/uppercase matching and compressed variants like `.fasta.gz`.
pub fn has_valid_extension<P: AsRef<Path>>(path: P, valid_extensions: &[&str], allow_compressed: bool) -> bool {
    let path = path.as_ref();

    // Get the outer extension (e.g., "gz" or "fasta")
    let ext = match path.extension().and_then(|e| e.to_str()) {
        Some(ext) => ext.to_ascii_lowercase(),
        None => return false,
    };

    if ext == "gz" && allow_compressed {
        // Strip .gz and inspect the underlying file extension
        let stem = match path.file_stem() {
            Some(s) => Path::new(s),
            None => return false,
        };

        stem.extension()
            .and_then(|e| e.to_str())
            .map(|e| valid_extensions.iter().any(|v| v.eq_ignore_ascii_case(e)))
            .unwrap_or(false)
    } else {
        valid_extensions.iter().any(|v| v.eq_ignore_ascii_case(&ext))
    }
}

/// caluclate gc percentage
#[inline]
pub fn get_gc_percentage(fasta_record: FastaRecord) -> f64 {
    let counts = fasta_record.num_atgc();
    let g_count = counts[2]; // G
    let c_count = counts[3]; // C
    let total_len = fasta_record.len();

    if total_len == 0 {
        0.0
    } else {
        ((g_count + c_count) as f64 / total_len as f64) * 100.0
    }
}

/// function to calculate N50 and N90
pub fn nx(lengths: &[usize], total: usize, fraction: f64) -> usize {
    let threshold = (total as f64 * fraction) as usize;

    let mut cumulative = 0;

    for &len in lengths {
        cumulative += len;

        if cumulative >= threshold {
            return len;
        }
    }

    0
}

/// function to format number 1000 to 1,000
pub fn num_to_str(mut n: u64) -> String {
    if n == 0 {
        return "0".to_string();
    }

    let mut result = String::new();
    let mut count = 0;

    while n > 0 {
        if count > 0 && count % 3 == 0 {
            result.push(',');
        }
        result.push((b'0' + (n % 10) as u8) as char);
        n /= 10;
        count += 1;
    }

    result.chars().rev().collect()
}

/// function to create a file
pub fn create_file<P: AsRef<Path>>(path: P) -> Result<File> {
    let p = path.as_ref();
    File::create(p).map_err(|e| AppError::FileIo {
        path: p.to_path_buf(),
        source: e,
    })
}

/// function to open a file
pub fn open_file<P: AsRef<Path>>(path: P) -> Result<File> {
    let p = path.as_ref();
    File::open(p).map_err(|e| AppError::FileIo {
        path: p.to_path_buf(),
        source: e,
    })
}

/// is compressed
/// Helper to determine if a path represents a compressed file based on its extension.
pub fn is_compressed<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_lowercase().as_str(), "gz" | "bgz" | "bz2" | "xz" | "zst"))
        .unwrap_or(false)
}

/// Formats a Duration into HH:MM:SS.mmm string format
pub fn format_time(duration: std::time::Duration) -> String {
    let total_millis = duration.as_millis();
    let millis = total_millis % 1000;
    let total_secs = total_millis / 1000;
    let seconds = total_secs % 60;
    let minutes = (total_secs / 60) % 60;
    let hours = total_secs / 3600;

    format!("{:02}:{:02}:{:02}.{:03}", hours, minutes, seconds, millis)
}