use crate::error::{AppError, Result};
use log::{debug, error, info};
use memmap2::Mmap;
use byteorder::{ByteOrder, LittleEndian};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use crate::utils::helpers::open_file;
use crate::data::kmer::Kmer;

const MAX_BYTE_COUNT: usize = 1 << 24; // 16 MB chunking page size

/// Minimal Kmer representation interface expected by the KMC lookup logic.
pub trait KmerRef {
    fn get_signature(&self, signature_len: usize) -> usize;
    fn get_prefix_fwd(&self) -> usize;
    fn get_suffix_fwd(&self) -> &[u8];
}

enum SuffixStorage {
    /// Thread-safe zero-copy memory-mapped chunks
    Mapped(Vec<Mmap>),
    /// Pre-loaded fully in-memory byte buffers
    InMemory(Vec<Vec<u8>>),
}

pub struct Kmc {
    kmc_prefix_file: PathBuf,
    kmc_suffix_file: PathBuf,
    kmer_length: u32,
    mode: u32,
    counter_size: usize,
    lut_prefix_length: usize,
    suffix_length: usize,
    signature_length: usize,
    min_count: u32,
    max_count: u32,
    total_kmers: u64,
    both_strands: bool,
    single_lut_size: u64,
    version: u32,
    record_size: usize,
    records_per_page: u64,
    lut_prefix_array_size: usize,
    prefix_array: Vec<u64>,
    signature_map: Vec<u32>,
    suffix_storage: SuffixStorage,
}

impl Kmc {
    pub fn new<P: AsRef<Path>>(kmc_db_name: P, in_memory: bool) -> Result<Self> {
        let prefix_path = PathBuf::from(format!("{}.kmc_pre", kmc_db_name.as_ref().display()));
        let suffix_path = PathBuf::from(format!("{}.kmc_suf", kmc_db_name.as_ref().display()));

        let mut kmc = Self::read_prefix_file(&prefix_path)?;
        kmc.kmc_prefix_file = prefix_path;
        kmc.kmc_suffix_file = suffix_path;

        if in_memory {
            kmc.preload_suffix_buffers()?;
        } else {
            kmc.read_suffix_buffers()?;
        }

        kmc.print_summary();
        Ok(kmc)
    }

    fn read_prefix_file(prefix_path: &Path) -> Result<Self> {
        info!("reading KMC prefix file {:?}", prefix_path.display());
        let file = open_file(prefix_path)?;
        let mmap = unsafe { Mmap::map(&file).map_err(AppError::Io)? };
        let file_size = mmap.len();

        if file_size < 12 {
            return Err(AppError::Generic("KMC prefix file too small".into()));
        }

        let header_offset = LittleEndian::read_u32(&mmap[file_size - 8..file_size - 4]) as usize;
        let header_start = file_size - header_offset - 8;
        let mut pos = header_start;

        let kmer_length = LittleEndian::read_u32(&mmap[pos..pos + 4]); pos += 4;
        let mode = LittleEndian::read_u32(&mmap[pos..pos + 4]); pos += 4;
        let counter_size = LittleEndian::read_u32(&mmap[pos..pos + 4]) as usize; pos += 4;
        let lut_prefix_length = LittleEndian::read_u32(&mmap[pos..pos + 4]) as usize; pos += 4;
        let suffix_length = (kmer_length as usize) - lut_prefix_length;
        let signature_length = LittleEndian::read_u32(&mmap[pos..pos + 4]) as usize; pos += 4;
        let min_count = LittleEndian::read_u32(&mmap[pos..pos + 4]); pos += 4;
        let max_count = LittleEndian::read_u32(&mmap[pos..pos + 4]); pos += 4;
        let total_kmers = LittleEndian::read_u64(&mmap[pos..pos + 8]); pos += 8;
        let both_strands = mmap[pos] == 0; pos += 1;

        pos += 3;  // skip uchar[3] padding
        pos += 24; // skip uint32[6] padding

        let version = LittleEndian::read_u32(&mmap[pos..pos + 4]);
        if version != 0x200 {
            error!("KMC version is not 0x200 (found 0x{:x})", version);
        }

        let signature_map_size = (1u64 << (2 * signature_length)) as usize + 1;
        let signature_map_start = file_size - header_offset - 8 - (signature_map_size * 4);

        let mut signature_map = Vec::with_capacity(signature_map_size);
        for i in 0..signature_map_size {
            let offset = signature_map_start + (i * 4);
            signature_map.push(LittleEndian::read_u32(&mmap[offset..offset + 4]));
        }

        let prefix_array_start = 4;
        let lut_prefix_array_size = 1 << (2 * lut_prefix_length);
        let single_lut_size = (lut_prefix_array_size * 8) as u64;
        let num_prefix_arrays = (signature_map_start - 8 - 4) / single_lut_size as usize;

        let total_prefix_elems = num_prefix_arrays * lut_prefix_array_size;
        let mut prefix_array = Vec::with_capacity(total_prefix_elems);
        for i in 0..total_prefix_elems {
            let offset = prefix_array_start + (i * 8);
            prefix_array.push(LittleEndian::read_u64(&mmap[offset..offset + 8]));
        }

        let record_size = counter_size + (suffix_length / 4);
        let records_per_page = (MAX_BYTE_COUNT / record_size) as u64;

        Ok(Self {
            kmc_prefix_file: PathBuf::new(),
            kmc_suffix_file: PathBuf::new(),
            kmer_length,
            mode,
            counter_size,
            lut_prefix_length,
            suffix_length,
            signature_length,
            min_count,
            max_count,
            total_kmers,
            both_strands,
            single_lut_size,
            version,
            record_size,
            records_per_page,
            lut_prefix_array_size,
            prefix_array,
            signature_map,
            suffix_storage: SuffixStorage::Mapped(Vec::new()),
        })
    }

    fn read_suffix_buffers(&mut self) -> Result<()> {
        info!("memory mapping KMC suffix file {:?}", self.kmc_suffix_file.display());
        let file = open_file(&self.kmc_suffix_file)?;

        let full_page_size = (MAX_BYTE_COUNT / self.record_size) * self.record_size;
        let total_bytes = self.total_kmers * self.record_size as u64;
        let number_of_pages = (total_bytes / full_page_size as u64
            + if total_bytes % full_page_size as u64 == 0 { 0 } else { 1 }) as usize;

        let mut maps = Vec::with_capacity(number_of_pages);
        for i in 0..number_of_pages {
            let page_size = if i == number_of_pages - 1 {
                let rem = (total_bytes % full_page_size as u64) as usize;
                if rem == 0 { full_page_size } else { rem }
            } else {
                full_page_size
            };

            let offset = (full_page_size as u64 * i as u64) + 4;
            let mmap = unsafe {
                memmap2::MmapOptions::new()
                    .offset(offset)
                    .len(page_size)
                    .map(&file)
                    .map_err(AppError::Io)?
            };
            maps.push(mmap);
        }

        self.suffix_storage = SuffixStorage::Mapped(maps);
        Ok(())
    }

    fn preload_suffix_buffers(&mut self) -> Result<()> {
        info!("loading KMC suffix file {:?} into memory", self.kmc_suffix_file.display());
        let mut file = File::open(&self.kmc_suffix_file).map_err(AppError::Io)?;

        let full_page_size = (MAX_BYTE_COUNT / self.record_size) * self.record_size;
        let total_bytes = self.total_kmers * self.record_size as u64;
        let number_of_pages = (total_bytes / full_page_size as u64
            + if total_bytes % full_page_size as u64 == 0 { 0 } else { 1 }) as usize;

        file.seek(SeekFrom::Start(4)).map_err(AppError::Io)?;

        let mut buffers = Vec::with_capacity(number_of_pages);
        for i in 0..number_of_pages {
            let page_size = if i == number_of_pages - 1 {
                let rem = (total_bytes % full_page_size as u64) as usize;
                if rem == 0 { full_page_size } else { rem }
            } else {
                full_page_size
            };

            let mut page_buf = vec![0u8; page_size];
            file.read_exact(&mut page_buf).map_err(AppError::Io)?;
            buffers.push(page_buf);
        }

        self.suffix_storage = SuffixStorage::InMemory(buffers);
        Ok(())
    }

    #[inline(always)]
    pub fn get_entry(&self, index: u64) -> &[u8] {
        // High-performance replacement for / and %
        let page = (index / self.records_per_page) as usize;
        let offset = (index % self.records_per_page) as usize;
        let record_offset = offset * self.record_size;

        match &self.suffix_storage {
            SuffixStorage::Mapped(maps) => {
                unsafe {
                    let page_slice = maps.get_unchecked(page);
                    page_slice.get_unchecked(record_offset..record_offset + self.record_size)
                }
            }
            SuffixStorage::InMemory(buffers) => {
                unsafe {
                    let page_slice = buffers.get_unchecked(page);
                    page_slice.get_unchecked(record_offset..record_offset + self.record_size)
                }
            }
        }
    }

    #[inline]
    pub fn get_suffix_from_entry<'a>(&self, entry: &'a [u8]) -> &'a [u8] {
        &entry[..self.suffix_length / 4]
    }

    #[inline]
    pub fn get_count_from_entry(&self, entry: &[u8]) -> u32 {
        let mut count = 0u32;
        let suffix_bytes = self.suffix_length / 4;
        for i in 0..self.counter_size {
            count |= (entry[suffix_bytes + i] as u32) << (i * 8);
        }
        count
    }

    /// Primary lookup function for k-mer frequencies.
    pub fn get_count<K: KmerRef>(&self, kmer: &K) -> u32 {
        let signature = kmer.get_signature(self.signature_length);
        let prefix = kmer.get_prefix_fwd();
        let suffix = kmer.get_suffix_fwd();

        let signature_offset = self.signature_map[signature] as usize;
        let start_idx = signature_offset * self.lut_prefix_array_size + prefix;

        if start_idx >= self.prefix_array.len() {
            return 0;
        }

        let start_pos = self.prefix_array[start_idx];
        let end_pos = if start_idx + 1 < self.prefix_array.len() {
            self.prefix_array[start_idx + 1]
        } else {
            self.total_kmers
        };

        // If start_pos >= end_pos, no suffixes exist for this prefix
        if start_pos >= end_pos {
            return 0;
        }

        let mut start = start_pos;
        let mut end = end_pos - 1;

        while start <= end {
            let mid = start + (end - start) / 2;
            let entry = self.get_entry(mid);
            let suffix_from_entry = self.get_suffix_from_entry(entry);

            match suffix.cmp(suffix_from_entry) {
                std::cmp::Ordering::Less => {
                    if mid == 0 {
                        break;
                    }
                    end = mid - 1;
                }
                std::cmp::Ordering::Greater => {
                    start = mid + 1;
                }
                std::cmp::Ordering::Equal => {
                    return self.get_count_from_entry(entry);
                }
            }
        }

        0
    }

    pub fn is_exist<K: KmerRef>(&self, kmer: &K) -> bool {
        self.get_count(kmer) != 0
    }

    pub fn convert_bytes_to_kmer(&self, suffix_bytes: &[u8], suffix_length: usize) -> String {
        let mut suffix = String::with_capacity(suffix_length);
        for i in 0..suffix_length {
            let base_index = (suffix_bytes[i / 4] >> (6 - 2 * (i % 4))) & 0x03;
            suffix.push(match base_index {
                0b00 => 'A',
                0b01 => 'C',
                0b10 => 'G',
                0b11 => 'T',
                _ => unreachable!(),
            });
        }
        suffix
    }

    pub fn print_summary(&self) {
        info!("==================== KMC INFO ====================");
        info!("{:<25}: {:?}", "KMC prefix file", self.kmc_prefix_file);
        info!("{:<25}: {:?}", "KMC suffix file", self.kmc_suffix_file);
        info!("{:<25}: {}", "Kmer length", self.kmer_length);
        info!("{:<25}: {}", "Mode", self.mode);
        info!("{:<25}: {}", "Counter size", self.counter_size);
        info!("{:<25}: {}", "LUT prefix length", self.lut_prefix_length);
        info!("{:<25}: {}", "Signature length", self.signature_length);
        info!("{:<25}: {}", "Min count", self.min_count);
        info!("{:<25}: {}", "Max count", self.max_count);
        info!("{:<25}: {}", "Total kmers", self.total_kmers);
        info!("{:<25}: {}", "Both strands", self.both_strands);
        info!("{:<25}: {}", "Signature map size", self.signature_map.len());
        info!("{:<25}: {}", "Single LUT size", self.single_lut_size);
        info!("{:<25}: {}", "Prefix buffer size", self.prefix_array.len());
        info!("{:<25}: {:x}", "Version", self.version);
        info!("==================================================");
    }

    pub fn dump_prefix_array<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let column_size = 1 << (2 * self.lut_prefix_length);
        info!("writing prefix array to file {:?}", path.as_ref().display());
        let file = File::create(path).map_err(AppError::Io)?;
        let mut writer = BufWriter::new(file);

        for (i, &val) in self.prefix_array.iter().enumerate() {
            write!(writer, "{}", val).map_err(AppError::Io)?;
            if i % column_size == column_size - 1 {
                writeln!(writer).map_err(AppError::Io)?;
            } else {
                write!(writer, "\t").map_err(AppError::Io)?;
            }
        }
        Ok(())
    }

    /// Optimized dump with page-aware batch streaming & fast 4-base unpacking
    pub fn dump<W: Write>(&self, mut writer: W) -> Result<()> {
        let total_prefix_length = self.prefix_array.len();

        let full_kmer_len = self.lut_prefix_length + self.suffix_length;
        // Allocate a single line buffer including '\t', count placeholder, and '\n'
        let mut line_buf = vec![b'A'; full_kmer_len];
        let mut itoa_buf = itoa::Buffer::new();

        const LOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];

        for i in 0..total_prefix_length {
            let start = self.prefix_array[i];
            let end = if i == total_prefix_length - 1 {
                self.total_kmers.saturating_sub(1)
            } else {
                self.prefix_array[i + 1].saturating_sub(1)
            };

            if start > end {
                continue;
            }

            // Decode the prefix directly into line_buf
            let prefix_idx = i % self.lut_prefix_array_size;
            let mut val = prefix_idx;
            for pos in (0..self.lut_prefix_length).rev() {
                line_buf[pos] = LOOKUP[val & 0x03];
                val >>= 2;
            }

            // Stream k-mers continuously
            for j in start..=end {
                // Inline entry retrieval (bypasses division when sequential)
                let entry = self.get_entry(j);
                let suffix_bytes = self.get_suffix_from_entry(entry);
                let count = self.get_count_from_entry(entry);

                // Fast 4-bases-per-byte unpacking loop
                let mut out_pos = self.lut_prefix_length;
                for &byte in suffix_bytes {
                    line_buf[out_pos]     = LOOKUP[((byte >> 6) & 0x03) as usize];
                    line_buf[out_pos + 1] = LOOKUP[((byte >> 4) & 0x03) as usize];
                    line_buf[out_pos + 2] = LOOKUP[((byte >> 2) & 0x03) as usize];
                    line_buf[out_pos + 3] = LOOKUP[(byte & 0x03) as usize];
                    out_pos += 4;
                }

                // Batch I/O writes into a single contiguous write_all block
                writer.write_all(&line_buf[..full_kmer_len]).map_err(AppError::Io)?;
                writer.write_all(b"\t").map_err(AppError::Io)?;
                writer.write_all(itoa_buf.format(count).as_bytes()).map_err(AppError::Io)?;
                writer.write_all(b"\n").map_err(AppError::Io)?;
            }
        }

        writer.flush().map_err(AppError::Io)?;
        Ok(())
    }

    #[inline]
    pub fn kmer_length(&self) -> usize { self.kmer_length as usize }

    #[inline]
    pub fn lut_prefix_length(&self) -> usize { self.lut_prefix_length }

    #[inline]
    pub fn signature_length(&self) -> usize { self.signature_length }

    /// Helper for querying a generic `Kmer` directly
    pub fn get_kmer_count(&self, kmer: &Kmer) -> u32 {
        let query_kmer = KmcQueryKmer::from_kmc(*kmer, self);
        self.get_count(&query_kmer)
    }
}



pub struct KmcQueryKmer<'a> {
    kmer: Kmer<'a>,
    lut_prefix_len: usize,
    signature_len: usize,
    packed_suffix: Vec<u8>,
}

impl<'a> KmcQueryKmer<'a> {
    pub fn from_kmc(kmer: Kmer<'a>, kmc: &Kmc) -> Self {
        let sequence = kmer.sequence();
        let lut_prefix_len = kmc.lut_prefix_length();
        let suffix_bytes = &sequence[lut_prefix_len..];

        // Pack ASCII suffix bytes into KMC's 2-bit format (4 bases per byte)
        let mut packed_suffix = Vec::with_capacity((suffix_bytes.len() + 3) / 4);
        for chunk in suffix_bytes.chunks(4) {
            let mut byte = 0u8;
            for (i, &b) in chunk.iter().enumerate() {
                let val = match b {
                    b'A' | b'a' => 0b00,
                    b'C' | b'c' => 0b01,
                    b'G' | b'g' => 0b10,
                    b'T' | b't' => 0b11,
                    _ => 0b00,
                };
                byte |= val << (6 - 2 * i);
            }
            packed_suffix.push(byte);
        }

        Self {
            kmer,
            lut_prefix_len,
            signature_len: kmc.signature_length(),
            packed_suffix,
        }
    }

    /// Exact KMC canonical signature value calculation for an S-mer
    #[inline]
    fn calc_signature_val(slice: &[u8]) -> usize {
        let mut sig = 0usize;
        for (i, &byte) in slice.iter().enumerate() {
            let val = match byte {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _ => 0b00,
            };
            sig |= (val as usize) << (2 * i);
        }
        sig
    }
}

impl<'a> KmerRef for KmcQueryKmer<'a> {
    #[inline]
    fn get_signature(&self, signature_len: usize) -> usize {
        if signature_len == 0 {
            return 0;
        }

        let seq = self.kmer.sequence();
        if seq.len() < signature_len {
            return 0;
        }

        let mut min_sig = usize::MAX;

        // Slide window of size signature_len across sequence
        for window in seq.windows(signature_len) {
            let sig = Self::calc_signature_val(window);

            // KMC ignores signature 0 (all A's) if signature_len > 1
            if sig < min_sig {
                min_sig = sig;
            }
        }

        if min_sig == usize::MAX { 0 } else { min_sig }
    }

    #[inline]
    fn get_prefix_fwd(&self) -> usize {
        let mut prefix = 0usize;
        let seq = self.kmer.sequence();
        let prefix_bytes = &seq[..self.lut_prefix_len];

        for &byte in prefix_bytes {
            let val = match byte {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _ => 0b00,
            };
            prefix = (prefix << 2) | val;
        }
        prefix
    }

    #[inline]
    fn get_suffix_fwd(&self) -> &[u8] {
        &self.packed_suffix
    }
}
