use crate::error::{AppError, Result};
use crate::utils::defaults::{DEFAULT_BUF_SIZE, DEFAULT_COLUMN_WIDTH};
use crate::utils::helpers::{create_file, num_to_str, open_file};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::{debug, error, info, trace};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Standard 4-bit BAM sequence decoder lookup table
const BAM_BASE_LOOKUP: [u8; 16] = [
    b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

/// Lightweight representation of a parsed BAM alignment record.
#[derive(Debug, Clone, Default)]
pub struct BamRecord {
    pub id: usize,
    pub name: String,
    pub flag: u16,
    pub ref_id: i32,
    pub pos: i32, // 0-based coordinate (-1 if unmapped)
    pub mapq: u8,
    pub next_ref_id: i32,
    pub next_pos: i32,
    pub tlen: i32,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,      // Raw Phred scores (0-255)
    pub nm: Option<i32>,       // Edit distance (NM auxiliary tag)
    pub as_score: Option<i32>, // Alignment score (AS auxiliary tag)
}

impl BamRecord {
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Calculate mean Phred quality score (raw values, no offset required).
    pub fn mean_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }
        let sum: usize = self.quality.iter().map(|&q| q as usize).sum();
        sum as f64 / self.quality.len() as f64
    }

    /// Count bases meeting a minimum Phred score threshold (e.g., Q20, Q30).
    pub fn count_q_bases(&self, min_q: u8) -> usize {
        self.quality.iter().filter(|&&q| q >= min_q).count()
    }

    /// Returns edit distance (mismatches/indels) if the NM tag is present.
    #[inline]
    pub fn nm(&self) -> Option<i32> {
        self.nm
    }

    /// Returns alignment score if the AS tag is present.
    #[inline]
    pub fn alignment_score(&self) -> Option<i32> {
        self.as_score
    }

    // -------------------------------------------------------------
    // SAM/BAM Bitwise FLAG Helpers
    // -------------------------------------------------------------

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.flag & 0x1 != 0
    }

    #[inline]
    pub fn is_proper_pair(&self) -> bool {
        self.flag & 0x2 != 0
    }

    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }

    #[inline]
    pub fn is_mapped(&self) -> bool {
        self.flag & 0x4 == 0
    }

    #[inline]
    pub fn mate_unmapped(&self) -> bool {
        self.flag & 0x8 != 0
    }

    #[inline]
    pub fn is_mate_unmapped(&self) -> bool {
        self.flag & 0x8 != 0
    }

    #[inline]
    pub fn is_mate_mapped(&self) -> bool {
        self.flag & 0x8 == 0
    }

    #[inline]
    pub fn is_first_in_pair(&self) -> bool {
        self.flag & 0x40 != 0
    }

    #[inline]
    pub fn is_second_in_pair(&self) -> bool {
        self.flag & 0x80 != 0
    }

    #[inline]
    pub fn is_secondary(&self) -> bool {
        self.flag & 0x100 != 0
    }

    #[inline]
    pub fn is_duplicate(&self) -> bool {
        self.flag & 0x400 != 0
    }

    #[inline]
    pub fn is_supplementary(&self) -> bool {
        self.flag & 0x800 != 0
    }
}

/// Represents a BAM reference sequence from the header block.
#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    pub name: String,
    pub length: i32,
}

/// Native BAM Header parsed from binary BAM stream.
#[derive(Debug, Clone, Default)]
pub struct BamHeader {
    pub text: String,
    pub references: Vec<ReferenceSequence>,
}

pub struct NativeBamReader<R> {
    reader: R,
    pub header: BamHeader,
    current_id: usize,
}

impl NativeBamReader<Box<dyn Read>> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();

        let file = open_file(path_ref)?;
        info!("reading native BAM file {:?}", path_ref.display());

        let file_buf = BufReader::with_capacity(DEFAULT_BUF_SIZE, file);
        let gz_decoder = MultiGzDecoder::new(file_buf);
        let reader: Box<dyn Read> = Box::new(BufReader::with_capacity(DEFAULT_BUF_SIZE, gz_decoder));

        Self::new(reader).map_err(Into::into)
    }
}

impl<R: Read> NativeBamReader<R> {
    pub fn new(mut reader: R) -> io::Result<Self> {
        // 1. Verify BAM Magic Bytes ("BAM\1")
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != b"BAM\x01" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM file magic header",
            ));
        }

        // 2. Read SAM Header Text
        let mut l_text_bytes = [0u8; 4];
        reader.read_exact(&mut l_text_bytes)?;
        let l_text = i32::from_le_bytes(l_text_bytes) as usize;

        let mut text_buf = vec![0u8; l_text];
        reader.read_exact(&mut text_buf)?;
        let text = String::from_utf8_lossy(&text_buf).into_owned();

        // 3. Read Reference Dictionary
        let mut n_ref_bytes = [0u8; 4];
        reader.read_exact(&mut n_ref_bytes)?;
        let n_ref = i32::from_le_bytes(n_ref_bytes) as usize;

        let mut references = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let mut l_name_bytes = [0u8; 4];
            reader.read_exact(&mut l_name_bytes)?;
            let l_name = i32::from_le_bytes(l_name_bytes) as usize;

            let mut name_buf = vec![0u8; l_name];
            reader.read_exact(&mut name_buf)?;
            let name = String::from_utf8_lossy(&name_buf[..l_name.saturating_sub(1)]).into_owned(); // Trim null terminator

            let mut l_ref_bytes = [0u8; 4];
            reader.read_exact(&mut l_ref_bytes)?;
            let length = i32::from_le_bytes(l_ref_bytes);

            references.push(ReferenceSequence { name, length });
        }

        let header = BamHeader { text, references };

        Ok(Self {
            reader,
            header,
            current_id: 0,
        })
    }

    /// Streaming BAM record deserializer without allocating helper structures.
    pub fn read_next(&mut self, record: &mut BamRecord) -> io::Result<bool> {
        let mut block_size_bytes = [0u8; 4];
        match self.reader.read_exact(&mut block_size_bytes) {
            Ok(_) => {}
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                debug!("reached end of BAM stream after {} records.", self.current_id);
                return Ok(false);
            }
            Err(e) => return Err(e),
        }

        let block_size = i32::from_le_bytes(block_size_bytes) as usize;
        let mut block = vec![0u8; block_size];
        self.reader.read_exact(&mut block)?;

        if block_size < 32 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "truncated BAM record header block",
            ));
        }

        // Parse fixed core alignment fields (32 bytes total)
        record.ref_id = i32::from_le_bytes(block[0..4].try_into().unwrap());
        record.pos = i32::from_le_bytes(block[4..8].try_into().unwrap());

        let l_read_name = block[8] as usize;
        record.mapq = block[9];

        let _bin = u16::from_le_bytes(block[10..12].try_into().unwrap());
        let n_cigar_op = u16::from_le_bytes(block[12..14].try_into().unwrap()) as usize;
        record.flag = u16::from_le_bytes(block[14..16].try_into().unwrap());

        let l_seq = i32::from_le_bytes(block[16..20].try_into().unwrap()) as usize;
        record.next_ref_id = i32::from_le_bytes(block[20..24].try_into().unwrap());
        record.next_pos = i32::from_le_bytes(block[24..28].try_into().unwrap());
        record.tlen = i32::from_le_bytes(block[28..32].try_into().unwrap());

        // Parse variable length fields
        let mut offset = 32;

        // 1. Read Name (null-terminated string)
        if offset + l_read_name > block_size {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid read name length"));
        }
        let name_bytes = &block[offset..offset + l_read_name.saturating_sub(1)];
        record.name = String::from_utf8_lossy(name_bytes).into_owned();
        offset += l_read_name;

        // 2. Skip CIGAR (n_cigar_op * 4 bytes)
        offset += n_cigar_op * 4;

        // 3. Decode 4-bit compressed nucleotide Sequence
        let seq_bytes_len = (l_seq + 1) / 2;
        if offset + seq_bytes_len > block_size {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid sequence length"));
        }

        record.sequence.clear();
        record.sequence.reserve(l_seq);
        for &byte in &block[offset..offset + seq_bytes_len] {
            let high = (byte >> 4) as usize;
            let low = (byte & 0x0F) as usize;
            record.sequence.push(BAM_BASE_LOOKUP[high]);
            if record.sequence.len() < l_seq {
                record.sequence.push(BAM_BASE_LOOKUP[low]);
            }
        }
        offset += seq_bytes_len;

        // 4. Quality Scores
        if offset + l_seq > block_size {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid quality length"));
        }
        record.quality.clear();
        record.quality.extend_from_slice(&block[offset..offset + l_seq]);
        offset += l_seq;

        // 5. Parse Auxiliary Tags (Extract NM, AS if available)
        record.nm = None;
        record.as_score = None;

        while offset + 3 <= block_size {
            let tag = [block[offset], block[offset + 1]];
            let val_type = block[offset + 2];
            offset += 3;

            let val = match val_type {
                b'c' => {
                    if offset + 1 > block_size { break; }
                    let v = block[offset] as i8 as i32;
                    offset += 1;
                    Some(v)
                }
                b'C' => {
                    if offset + 1 > block_size { break; }
                    let v = block[offset] as u8 as i32;
                    offset += 1;
                    Some(v)
                }
                b's' => {
                    if offset + 2 > block_size { break; }
                    let v = i16::from_le_bytes(block[offset..offset + 2].try_into().unwrap()) as i32;
                    offset += 2;
                    Some(v)
                }
                b'S' => {
                    if offset + 2 > block_size { break; }
                    let v = u16::from_le_bytes(block[offset..offset + 2].try_into().unwrap()) as i32;
                    offset += 2;
                    Some(v)
                }
                b'i' => {
                    if offset + 4 > block_size { break; }
                    let v = i32::from_le_bytes(block[offset..offset + 4].try_into().unwrap());
                    offset += 4;
                    Some(v)
                }
                b'I' => {
                    if offset + 4 > block_size { break; }
                    let v = u32::from_le_bytes(block[offset..offset + 4].try_into().unwrap()) as i32;
                    offset += 4;
                    Some(v)
                }
                b'f' => {
                    if offset + 4 > block_size { break; }
                    offset += 4;
                    None
                }
                b'A' => {
                    if offset + 1 > block_size { break; }
                    offset += 1;
                    None
                }
                b'Z' | b'H' => {
                    while offset < block_size && block[offset] != 0 {
                        offset += 1;
                    }
                    if offset < block_size {
                        offset += 1; // skip null byte
                    }
                    None
                }
                b'B' => {
                    if offset + 5 > block_size { break; }
                    let elem_type = block[offset];
                    let count = u32::from_le_bytes(block[offset + 1..offset + 5].try_into().unwrap()) as usize;
                    offset += 5;
                    let elem_size = match elem_type {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        b'i' | b'I' | b'f' => 4,
                        _ => 1,
                    };
                    offset += count * elem_size;
                    None
                }
                _ => break,
            };

            if tag == *b"NM" {
                record.nm = val;
            } else if tag == *b"AS" {
                record.as_score = val;
            }
        }

        record.id = self.current_id;
        self.current_id += 1;

        trace!("parsed BAM record ID {}: name='{}', len={} bp", record.id, record.name, record.len());
        Ok(true)
    }
}

impl<R: Read> Iterator for NativeBamReader<R> {
    type Item = io::Result<BamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = BamRecord::default();
        match self.read_next(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct BamFlagStats {
    pub total: u64,
    pub primary: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub duplicates: u64,
    pub primary_duplicates: u64,
    pub mapped: u64,
    pub primary_mapped: u64,
    pub paired_in_seq: u64,
    pub read1: u64,
    pub read2: u64,
    pub properly_paired: u64,
    pub mate_mapped: u64,
    pub singletons: u64,
    pub diff_chr: u64,
    pub diff_chr_mapq5: u64,
}

impl BamFlagStats {
    pub fn update(&mut self, rec: &BamRecord) {
        self.total += 1;

        if rec.is_secondary() {
            self.secondary += 1;
        } else if rec.is_supplementary() {
            self.supplementary += 1;
        } else {
            self.primary += 1;
        }

        if rec.is_duplicate() {
            self.duplicates += 1;
            if !rec.is_secondary() && !rec.is_supplementary() {
                self.primary_duplicates += 1;
            }
        }

        if rec.is_mapped() {
            self.mapped += 1;
            if !rec.is_secondary() && !rec.is_supplementary() {
                self.primary_mapped += 1;
            }
        }

        if rec.is_paired() {
            self.paired_in_seq += 1;
            if rec.is_first_in_pair() {
                self.read1 += 1;
            }
            if rec.is_second_in_pair() {
                self.read2 += 1;
            }
            if rec.is_proper_pair() && rec.is_mapped() {
                self.properly_paired += 1;
            }
            if rec.is_mapped() && rec.is_mate_mapped() {
                self.mate_mapped += 1;
                if rec.ref_id != rec.next_ref_id && rec.ref_id >= 0 && rec.next_ref_id >= 0 {
                    self.diff_chr += 1;
                    if rec.mapq >= 5 {
                        self.diff_chr_mapq5 += 1;
                    }
                }
            } else if rec.is_mapped() && rec.is_mate_unmapped() {
                self.singletons += 1;
            }
        }
    }

    pub fn log_summary(&self) {
        let pct = |val: u64, base: u64| -> String {
            if base > 0 {
                format!("{:.2}%", (val as f64 / base as f64) * 100.0)
            } else {
                "0.00%".to_string()
            }
        };

        let metrics = vec![
            ("Total reads", num_to_str(self.total)),
            ("Primary reads", num_to_str(self.primary)),
            ("Secondary alignments", num_to_str(self.secondary)),
            ("Supplementary alignments", num_to_str(self.supplementary)),
            ("Duplicate reads", num_to_str(self.duplicates)),
            (
                "Mapped reads",
                format!("{} ({})", num_to_str(self.mapped), pct(self.mapped, self.total)),
            ),
            ("Paired in sequencing", num_to_str(self.paired_in_seq)),
            ("Read 1 count", num_to_str(self.read1)),
            ("Read 2 count", num_to_str(self.read2)),
            (
                "Properly paired",
                format!(
                    "{} ({})",
                    num_to_str(self.properly_paired),
                    pct(self.properly_paired, self.paired_in_seq)
                ),
            ),
            (
                "Both mates mapped",
                format!(
                    "{} ({})",
                    num_to_str(self.mate_mapped),
                    pct(self.mate_mapped, self.paired_in_seq)
                ),
            ),
            (
                "Singletons",
                format!(
                    "{} ({})",
                    num_to_str(self.singletons),
                    pct(self.singletons, self.paired_in_seq)
                ),
            ),
            (
                "Mate mapped to diff chr",
                format!(
                    "{} ({})",
                    num_to_str(self.diff_chr),
                    pct(self.diff_chr, self.paired_in_seq)
                ),
            ),
            (
                "Diff chr (MAPQ >= 5)",
                format!(
                    "{} ({})",
                    num_to_str(self.diff_chr_mapq5),
                    pct(self.diff_chr_mapq5, self.paired_in_seq)
                ),
            ),
        ];

        info!("=================== FLAGSTAT SUMMARY ===================");
        for (label, val) in metrics {
            info!("{:<width$} : {}", label, val, width = DEFAULT_COLUMN_WIDTH);
        }
        info!("========================================================");
    }
}