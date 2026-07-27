use crate::error::{AppError, Result};
use crate::utils::defaults::COMPLEMENT;
use crate::utils::defaults::{DEFAULT_BUF_SIZE, VALID_FASTQ_EXTENSIONS};
use crate::utils::helpers::create_file;
use crate::utils::helpers::{has_valid_extension, open_file, trim_ascii_whitespace};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::{debug, info, trace};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Represents a FASTQ record using raw byte vectors.
#[derive(Debug, Clone, Default)]
pub struct FastqRecord {
    pub id: usize,
    pub name: String,
    pub description: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl FastqRecord {
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Calculate mean Phred quality score (assuming Phred+33 offset).
    pub fn mean_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }
        let sum: usize = self.quality.iter().map(|&q| (q.saturating_sub(33)) as usize).sum();
        sum as f64 / self.quality.len() as f64
    }

    /// Count bases meeting a minimum Phred score threshold (e.g., Q20, Q30).
    pub fn count_q_bases(&self, min_q: u8) -> usize {
        let threshold = min_q + 33;
        self.quality.iter().filter(|&&q| q >= threshold).count()
    }

    /// In-place reverse complement of sequence AND reverse of quality scores.
    pub fn revcomp_mut(&mut self) {
        let len = self.sequence.len();
        for i in 0..len / 2 {
            let j = len - 1 - i;

            let a = COMPLEMENT[self.sequence[i] as usize];
            let b = COMPLEMENT[self.sequence[j] as usize];
            self.sequence[i] = b;
            self.sequence[j] = a;

            self.quality.swap(i, j);
        }
        if len % 2 != 0 {
            self.sequence[len / 2] = COMPLEMENT[self.sequence[len / 2] as usize];
        }
    }
}

pub struct FastqReader<R> {
    reader: R,
    line_buf: Vec<u8>,
    current_id: usize,
}

impl FastqReader<Box<dyn BufRead>> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();

        if !has_valid_extension(path_ref, VALID_FASTQ_EXTENSIONS, true) {
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

        info!("reading fastq file {:?} (gzipped: {})", path_ref.display(), is_gzipped);

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

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: Vec::with_capacity(512),
            current_id: 0,
        }
    }

    /// Zero-allocation record streaming into a mutable FASTQ target struct.
    pub fn read_next(&mut self, record: &mut FastqRecord) -> io::Result<bool> {
        // Line 1: Header (@name desc)
        self.line_buf.clear();
        let bytes_read = self.reader.read_until(b'\n', &mut self.line_buf)?;
        if bytes_read == 0 {
            debug!("reached end of FASTQ stream after {} records.", self.current_id);
            return Ok(false);
        }

        let trimmed_header = trim_ascii_whitespace(&self.line_buf);
        if !trimmed_header.starts_with(b"@") {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("corrupted FASTQ record {}: missing '@' prefix", self.current_id),
            ));
        }

        let mut parts = trimmed_header[1..].splitn(2, |&b| b == b' ' || b == b'\t');
        record.id = self.current_id;
        record.name = parts
            .next()
            .map(|b| String::from_utf8_lossy(b).into_owned())
            .unwrap_or_default();
        record.description = parts
            .next()
            .map(|b| String::from_utf8_lossy(trim_ascii_whitespace(b)).into_owned())
            .unwrap_or_default();

        // Line 2: Sequence
        self.line_buf.clear();
        if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "incomplete FASTQ record (missing sequence)",
            ));
        }
        record.sequence.clear();
        record.sequence.extend_from_slice(trim_ascii_whitespace(&self.line_buf));

        // Line 3: Separator (+)
        self.line_buf.clear();
        if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "incomplete FASTQ record (missing '+')",
            ));
        }

        // Line 4: Quality Scores
        self.line_buf.clear();
        if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "incomplete FASTQ record (missing quality)",
            ));
        }
        record.quality.clear();
        record.quality.extend_from_slice(trim_ascii_whitespace(&self.line_buf));

        if record.sequence.len() != record.quality.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FASTQ length mismatch for record '{}': sequence length {} != quality length {}",
                    record.name,
                    record.sequence.len(),
                    record.quality.len()
                ),
            ));
        }

        trace!("parsed FASTQ record ID {}: name='{}', len={} bp", record.id, record.name, record.len());

        self.current_id += 1;
        Ok(true)
    }
}

// Enable standard `for record in reader` loops
impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = FastqRecord::default();
        match self.read_next(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct PairedFastqReader<R1, R2> {
    r1_reader: FastqReader<R1>,
    r2_reader: FastqReader<R2>,
}

impl PairedFastqReader<Box<dyn BufRead>, Box<dyn BufRead>> {
    pub fn from_paths<P: AsRef<Path>>(r1_path: P, r2_path: P) -> Result<Self> {
        let r1_reader = FastqReader::from_path(r1_path)?;
        let r2_reader = FastqReader::from_path(r2_path)?;
        Ok(Self { r1_reader, r2_reader })
    }
}

impl<R1: BufRead, R2: BufRead> PairedFastqReader<R1, R2> {
    /// Zero-allocation paired streaming into pre-allocated `FastqRecord` target structs.
    pub fn read_next(&mut self, rec1: &mut FastqRecord, rec2: &mut FastqRecord) -> io::Result<bool> {
        let has_r1 = self.r1_reader.read_next(rec1)?;
        let has_r2 = self.r2_reader.read_next(rec2)?;

        match (has_r1, has_r2) {
            (true, true) => Ok(true),
            (false, false) => Ok(false),
            (true, false) => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "R1 file has more records than R2 file (unbalanced paired-end files)",
            )),
            (false, true) => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "R2 file has more records than R1 file (unbalanced paired-end files)",
            )),
        }
    }
}

// Enable standard `for pair in paired_reader` loops
impl<R1: BufRead, R2: BufRead> Iterator for PairedFastqReader<R1, R2> {
    type Item = io::Result<(FastqRecord, FastqRecord)>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut rec1 = FastqRecord::default();
        let mut rec2 = FastqRecord::default();
        match self.read_next(&mut rec1, &mut rec2) {
            Ok(true) => Some(Ok((rec1, rec2))),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct FastqWriter<W: Write> {
    writer: W,
}

impl FastqWriter<Box<dyn Write>> {
    pub fn create_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let file = create_file(path_ref)?;

        let is_gzipped = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);

        info!("created output FASTQ file {:?} (gzipped: {})", path_ref.display(), is_gzipped);

        let inner_buf = BufWriter::with_capacity(DEFAULT_BUF_SIZE, file);

        let writer: Box<dyn Write> = if is_gzipped {
            let gz_encoder = GzEncoder::new(inner_buf, Compression::fast());
            // Double-buffer wrapping GzEncoder with BufWriter for high throughput
            Box::new(BufWriter::with_capacity(DEFAULT_BUF_SIZE, gz_encoder))
        } else {
            Box::new(inner_buf)
        };

        Ok(Self::new(writer))
    }
}

impl<W: Write> FastqWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    pub fn write_record(&mut self, record: &FastqRecord) -> io::Result<()> {
        self.writer.write_all(b"@")?;
        self.writer.write_all(record.name.as_bytes())?;
        if !record.description.is_empty() {
            self.writer.write_all(b" ")?;
            self.writer.write_all(record.description.as_bytes())?;
        }
        self.writer.write_all(b"\n")?;
        self.writer.write_all(&record.sequence)?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(&record.quality)?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}