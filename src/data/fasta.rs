use crate::error::{AppError, Result};
use crate::utils::helpers::{create_file, has_valid_extension, open_file, trim_ascii_whitespace};
use crate::utils::defaults::{DEFAULT_BUF_SIZE, VALID_FASTA_EXTENSIONS, COMPLEMENT};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::{debug, info, trace}; // Added logging imports
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write, BufWriter};
use std::path::Path;
use crate::data::kmer::KmerIterator;

/// Represents a FASTA record using raw byte vectors to avoid UTF-8 validation overhead.
#[derive(Debug, Clone, Default)]
pub struct FastaRecord {
    pub id: usize,
    pub name: String,
    pub description: String,
    pub sequence: Vec<u8>,
}

impl FastaRecord {
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Optimized count using a fixed 256-element lookup array (zero branching).
    /// Returns counts for [A, T, G, C, N] (case-insensitive).
    pub fn num_atgc(&self) -> [usize; 5] {
        static LOOKUP: [u8; 256] = {
            let mut table = [255u8; 256];
            table[b'A' as usize] = 0; table[b'a' as usize] = 0;
            table[b'T' as usize] = 1; table[b't' as usize] = 1;
            table[b'G' as usize] = 2; table[b'g' as usize] = 2;
            table[b'C' as usize] = 3; table[b'c' as usize] = 3;
            table[b'N' as usize] = 4; table[b'n' as usize] = 4;
            table
        };

        let mut counts = [0usize; 5];
        for &byte in &self.sequence {
            let idx = LOOKUP[byte as usize];
            if idx != 255 {
                counts[idx as usize] += 1;
            }
        }
        counts
    }

    #[inline]
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }

    /// Calculate the GC percentage of the sequence.
    #[inline]
    pub fn gc_percentage(&self) -> f64 {
        crate::utils::helpers::get_gc_percentage(self.clone())
    }

    /// Get number of N's
    #[inline]
    pub fn num_n(&self) -> usize {
        let counts = self.num_atgc();
        counts[4] // N count
    }

    #[inline]
    pub fn revcomp(&self) -> Vec<u8> {
        let mut rev = Vec::with_capacity(self.sequence.len());
        for &byte in self.sequence.iter().rev() {
            // Unsafe lookup bypasses unnecessary bounds checks
            rev.push(unsafe { *COMPLEMENT.get_unchecked(byte as usize) });
        }
        rev
    }

    /// In-place reverse complement (Zero allocations!)
    pub fn revcomp_mut(&mut self) {
        let len = self.sequence.len();
        for i in 0..len / 2 {
            let j = len - 1 - i;
            let a = COMPLEMENT[self.sequence[i] as usize];
            let b = COMPLEMENT[self.sequence[j] as usize];
            self.sequence[i] = b;
            self.sequence[j] = a;
        }
        if len % 2 != 0 {
            self.sequence[len / 2] = COMPLEMENT[self.sequence[len / 2] as usize];
        }
    }

    /// Iterates over all k-mers of length `k`.
    #[inline]
    pub fn kmers(&self, k: usize) -> KmerIterator<'_> {
        KmerIterator::new(&self.sequence, k, false)
    }

    /// Iterates ONLY over canonical k-mers (lexicographically smaller strand).
    #[inline]
    pub fn canonical_kmers(&self, k: usize) -> KmerIterator<'_> {
        KmerIterator::new(&self.sequence, k, true)
    }
}

pub struct FastaReader<R> {
    reader: R,
    buffer: Vec<u8>,
    line_buf: Vec<u8>,
    current_id: usize,
    has_next_header: bool,
}

impl FastaReader<Box<dyn BufRead>> {
    /// Opens a file with standard or gzipped FASTA contents.
    /// Returns a trait object `Box<dyn BufRead>` to support both transparently.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();

        if !has_valid_extension(path_ref, VALID_FASTA_EXTENSIONS, true) {
            let ext = path_ref
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("none")
                .to_string();

            return Err(AppError::InvalidExtension {
                ext,
                path: path_ref.to_path_buf(),
            });
        }

        let file = open_file(path_ref)?;

        let is_gzipped = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);

        info!(
            "reading fasta file {:?} (gzipped: {})",
            path_ref.display(),
            is_gzipped
        );

        let reader: Box<dyn BufRead> = if is_gzipped {
            let file_buf = BufReader::with_capacity(DEFAULT_BUF_SIZE, file);
            let gz_decoder = MultiGzDecoder::new(file_buf);
            Box::new(BufReader::with_capacity(DEFAULT_BUF_SIZE, gz_decoder))
        } else {
            Box::new(BufReader::with_capacity(DEFAULT_BUF_SIZE, file))
        };

        Ok(Self::new(reader))
    }
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: Vec::with_capacity(1024),
            line_buf: Vec::with_capacity(256),
            current_id: 0,
            has_next_header: false,
        }
    }

    /// Zero-allocation record streaming into a mutable target struct.
    pub fn read_next(&mut self, record: &mut FastaRecord) -> io::Result<bool> {
        if !self.has_next_header {
            loop {
                self.line_buf.clear();
                let bytes_read = self.reader.read_until(b'\n', &mut self.line_buf)?;
                if bytes_read == 0 {
                    debug!("reached end of input stream after {} records.", self.current_id);
                    return Ok(false);
                }

                let trimmed = trim_ascii_whitespace(&self.line_buf);
                if trimmed.starts_with(b">") {
                    break;
                }
            }
        }

        let header = &self.line_buf[1..];
        let trimmed_header = trim_ascii_whitespace(header);

        let mut parts = trimmed_header.splitn(2, |&b| b == b' ' || b == b'\t');

        record.id = self.current_id;
        record.name = parts
            .next()
            .map(|b| String::from_utf8_lossy(b).into_owned())
            .unwrap_or_default();

        record.description = parts
            .next()
            .map(|b| String::from_utf8_lossy(trim_ascii_whitespace(b)).into_owned())
            .unwrap_or_default();

        record.sequence.clear();
        self.has_next_header = false;

        loop {
            self.line_buf.clear();
            let bytes_read = self.reader.read_until(b'\n', &mut self.line_buf)?;
            if bytes_read == 0 {
                break;
            }

            let trimmed = trim_ascii_whitespace(&self.line_buf);
            if trimmed.is_empty() {
                continue;
            }

            if trimmed.starts_with(b">") {
                self.has_next_header = true;
                break;
            }

            record.sequence.extend_from_slice(trimmed);
        }

        trace!(
            "parsed record ID {}: name='{}', len={} bp",
            record.id,
            record.name,
            record.len()
        );

        self.current_id += 1;
        Ok(true)
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = FastaRecord::default();
        match self.read_next(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct FastaWriter<W: Write> {
    writer: W,
}

impl FastaWriter<Box<dyn Write>> {
    pub fn create_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let file = create_file(path_ref)?;

        let is_gzipped = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);

        info!(
        "created output FASTA file {:?} (gzipped: {})",
        path_ref.display(),
        is_gzipped
    );

        let inner_buf = BufWriter::with_capacity(DEFAULT_BUF_SIZE, file);

        let writer: Box<dyn Write> = if is_gzipped {
            let gz_encoder = GzEncoder::new(inner_buf, Compression::fast());
            // CRITICAL FIX: Put a BufWriter ABOVE GzEncoder to batch small writes!
            Box::new(BufWriter::with_capacity(DEFAULT_BUF_SIZE, gz_encoder))
        } else {
            Box::new(inner_buf)
        };

        Ok(Self::new(writer))
    }
}

impl<W: Write> FastaWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    pub fn write_record(&mut self, record: &FastaRecord, fold: Option<usize>) -> io::Result<()> {
        self.write_raw(&record.name, Some(&record.description), &record.sequence, fold)
    }

    pub fn write_raw(
        &mut self,
        name: &str,
        description: Option<&str>,
        sequence: &[u8],
        fold: Option<usize>,
    ) -> io::Result<()> {
        trace!("writing record '{}", name);

        self.writer.write_all(b">")?;
        self.writer.write_all(name.as_bytes())?;

        if let Some(desc) = description {
            if !desc.is_empty() {
                self.writer.write_all(b" ")?;
                self.writer.write_all(desc.as_bytes())?;
            }
        }
        self.writer.write_all(b"\n")?;

        match fold {
            Some(width) if width > 0 => {
                for chunk in sequence.chunks(width) {
                    self.writer.write_all(chunk)?;
                    self.writer.write_all(b"\n")?;
                }
            }
            _ => {
                self.writer.write_all(sequence)?;
                self.writer.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    pub fn flush(&mut self) -> io::Result<()> {
        debug!("flushing output writer stream buffer");
        self.writer.flush()
    }
}