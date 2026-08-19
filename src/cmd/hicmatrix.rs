use crate::data::bam::{BamFlagStats, BamRecord, NativeBamReader};
use crate::error::{AppError, Result};
use crate::utils::defaults::DEFAULT_COLUMN_WIDTH;
use crate::utils::helpers::{create_file, num_to_str};
use crossbeam_channel::{bounded, Receiver, Sender};
use log::{info, warn};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::thread;

const CHUNK_SIZE: usize = 20_000;

/// Configurable options for Hi-C quality control read pair filtering.
#[derive(Debug, Clone, Copy)]
pub struct FilterOptions {
    pub min_mapq: u8,
    pub max_nm: Option<i32>,
    pub min_as: Option<i32>,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            min_mapq: 30,
            max_nm: Some(5),
            min_as: Some(100),
        }
    }
}

/// Represents a single genomic paired interaction record extracted from consecutive R1 & R2 BAM reads.
#[derive(Debug, Clone, Default)]
pub struct HiCContact {
    pub chrom1: String,
    pub pos1: usize,
    pub chrom2: String,
    pub pos2: usize,
}

/// Sparse matrix entry key representing a binned interaction pair.
#[derive(Hash, Eq, PartialEq, Clone, Debug)]
pub struct ContactBinKey {
    pub chrom1: String,
    pub bin1: usize, // 0-based bin start coordinate
    pub chrom2: String,
    pub bin2: usize,
}

impl ContactBinKey {
    /// Constructs a canonical bin pair key where (chrom1, bin1) <= (chrom2, bin2)
    /// to maintain upper-triangular sparse matrix symmetry.
    #[inline]
    pub fn new(mut c1: String, pos1: usize, mut c2: String, pos2: usize, bin_size: usize) -> Self {
        let mut b1 = (pos1 / bin_size) * bin_size;
        let mut b2 = (pos2 / bin_size) * bin_size;

        if c1 > c2 || (c1 == c2 && b1 > b2) {
            std::mem::swap(&mut c1, &mut c2);
            std::mem::swap(&mut b1, &mut b2);
        }

        Self {
            chrom1: c1,
            bin1: b1,
            chrom2: c2,
            bin2: b2,
        }
    }
}

/// Local worker state containing intermediate bin matrix counts and statistics.
#[derive(Default)]
struct MatrixChunkStats {
    total_pairs: usize,
    intra_pairs: usize,
    inter_pairs: usize,
    bin_counts: HashMap<ContactBinKey, u64>,
}

impl MatrixChunkStats {
    fn new() -> Self {
        Self {
            bin_counts: HashMap::with_capacity(10_000),
            ..Default::default()
        }
    }

    #[inline]
    fn merge(mut self, other: Self) -> Self {
        self.total_pairs += other.total_pairs;
        self.intra_pairs += other.intra_pairs;
        self.inter_pairs += other.inter_pairs;

        for (key, count) in other.bin_counts {
            *self.bin_counts.entry(key).or_insert(0) += count;
        }

        self
    }
}

/// Runs parallel Hi-C contact matrix binning using default filter options (min MAPQ = 30).
pub fn run_default(input_path: &str, bin_size: usize, nthreads: usize) -> Result<()> {
    run(input_path, bin_size, nthreads, FilterOptions::default())
}

/// Runs parallel Hi-C contact matrix binning by reading BAM records sequentially in pairs.
pub fn run(
    input_path: &str,
    bin_size: usize,
    nthreads: usize,
    filters: FilterOptions,
) -> Result<()> {
    info!(
        "starting Hi-C consecutive-pair contact matrix binning (bin size: {} bp, file: {}, threads: {}, min MAPQ: {}, max NM: {:?}, min AS: {:?})",
        num_to_str(bin_size as u64),
        input_path,
        nthreads,
        filters.min_mapq,
        filters.max_nm,
        filters.min_as
    );

    let worker_threads = if nthreads > 1 { nthreads - 1 } else { 1 };

    let (tx_chunks, rx_chunks): (Sender<Vec<HiCContact>>, Receiver<Vec<HiCContact>>) = bounded(100);
    let (tx_stats, rx_stats): (Sender<MatrixChunkStats>, Receiver<MatrixChunkStats>) = bounded(100);

    // -------------------------------------------------------------
    // 1. PRODUCER THREAD (Pair-wise BAM Streaming Reader)
    // -------------------------------------------------------------
    let input_clone = input_path.to_string();
    let reader_handle = thread::spawn(move || -> Result<BamFlagStats> {
        let mut bam_reader = NativeBamReader::from_path(&input_clone)?;
        let mut flag_stats = BamFlagStats::default();

        let ref_names: Vec<String> = bam_reader
            .header
            .references
            .iter()
            .map(|r| r.name.clone())
            .collect();

        let mut current_rec = BamRecord::default();
        let mut pending_rec: Option<BamRecord> = None;
        let mut chunk_buffer = Vec::with_capacity(CHUNK_SIZE);

        while bam_reader.read_next(&mut current_rec)? {
            flag_stats.update(&current_rec);

            // Ignore non-primary alignments immediately to prevent stream desynchronization
            if current_rec.is_secondary() || current_rec.is_supplementary() {
                continue;
            }

            if let Some(prev) = pending_rec.take() {
                if prev.name == current_rec.name {
                    let rec1 = prev;
                    let rec2 = current_rec;
                    current_rec = BamRecord::default();

                    // Quality control filtering on the primary pair
                    if rec1.is_unmapped()
                        || rec2.is_unmapped()
                        || rec1.ref_id < 0
                        || rec2.ref_id < 0
                        || rec1.is_duplicate()
                        || rec2.is_duplicate()
                        || rec1.mapq < filters.min_mapq
                        || rec2.mapq < filters.min_mapq
                    {
                        continue;
                    }

                    // Filter on maximum allowable edit distance (NM tag)
                    if let Some(max_nm) = filters.max_nm {
                        if rec1.nm.map_or(false, |nm| nm > max_nm)
                            || rec2.nm.map_or(false, |nm| nm > max_nm)
                        {
                            continue;
                        }
                    }

                    // Filter on minimum allowable alignment score (AS tag)
                    if let Some(min_as) = filters.min_as {
                        if rec1.as_score.map_or(false, |score| score < min_as)
                            || rec2.as_score.map_or(false, |score| score < min_as)
                        {
                            continue;
                        }
                    }

                    let r1_idx = rec1.ref_id as usize;
                    let r2_idx = rec2.ref_id as usize;

                    if r1_idx < ref_names.len() && r2_idx < ref_names.len() {
                        chunk_buffer.push(HiCContact {
                            chrom1: ref_names[r1_idx].clone(),
                            pos1: rec1.pos as usize,
                            chrom2: ref_names[r2_idx].clone(),
                            pos2: rec2.pos as usize,
                        });
                    }

                    if chunk_buffer.len() >= CHUNK_SIZE {
                        let full_chunk = std::mem::take(&mut chunk_buffer);
                        if tx_chunks.send(full_chunk).is_err() {
                            break;
                        }
                        chunk_buffer.reserve(CHUNK_SIZE);
                    }
                } else {
                    // Name mismatch: prev was a singleton or missing mate. Keep current_rec for next pair.
                    pending_rec = Some(current_rec);
                    current_rec = BamRecord::default();
                }
            } else {
                pending_rec = Some(current_rec);
                current_rec = BamRecord::default();
            }
        }

        if !chunk_buffer.is_empty() {
            let _ = tx_chunks.send(chunk_buffer);
        }

        Ok(flag_stats)
    });

    // -------------------------------------------------------------
    // 2. CONSUMER THREADS (Parallel Binning & Local Accumulation)
    // -------------------------------------------------------------
    let mut worker_handles = Vec::with_capacity(worker_threads);

    for _ in 0..worker_threads {
        let rx = rx_chunks.clone();
        let tx = tx_stats.clone();

        let handle = thread::spawn(move || {
            let mut local_stats = MatrixChunkStats::new();

            while let Ok(chunk) = rx.recv() {
                for contact in chunk {
                    process_single_contact(&contact, bin_size, &mut local_stats);
                }
            }

            let _ = tx.send(local_stats);
        });

        worker_handles.push(handle);
    }

    drop(rx_chunks);
    drop(tx_stats);

    // -------------------------------------------------------------
    // 3. AGGREGATOR (Main Thread - Matrix Merge & Thread Joining)
    // -------------------------------------------------------------
    let mut combined_stats = MatrixChunkStats::new();
    while let Ok(stats) = rx_stats.recv() {
        combined_stats = combined_stats.merge(stats);
    }

    let flag_stats = reader_handle.join().unwrap()?;
    for handle in worker_handles {
        handle.join().unwrap();
    }

    // -------------------------------------------------------------
    // 4. SUMMARY LOGGING & TSV MATRIX WRITING
    // -------------------------------------------------------------
    flag_stats.log_summary();

    let summary = CalculatedMatrixSummary::from_stats(input_path, bin_size, &filters, &combined_stats);
    summary.log_summary();
    summary.write_tsv(&combined_stats.bin_counts)?;

    Ok(())
}

#[inline]
fn process_single_contact(contact: &HiCContact, bin_size: usize, stats: &mut MatrixChunkStats) {
    stats.total_pairs += 1;

    if contact.chrom1 == contact.chrom2 {
        stats.intra_pairs += 1;
    } else {
        stats.inter_pairs += 1;
    }

    let bin_key = ContactBinKey::new(
        contact.chrom1.clone(),
        contact.pos1,
        contact.chrom2.clone(),
        contact.pos2,
        bin_size,
    );

    *stats.bin_counts.entry(bin_key).or_insert(0) += 1;
}

struct CalculatedMatrixSummary<'a> {
    fname: &'a str,
    bin_size: usize,
    filters: &'a FilterOptions,
    total_pairs: usize,
    intra_pairs: usize,
    inter_pairs: usize,
    populated_bins: usize,
    intra_pct: f64,
    inter_pct: f64,
}

impl<'a> CalculatedMatrixSummary<'a> {
    fn from_stats(
        fname: &'a str,
        bin_size: usize,
        filters: &'a FilterOptions,
        stats: &MatrixChunkStats,
    ) -> Self {
        let total_pairs = stats.total_pairs;
        let (intra_pct, inter_pct) = if total_pairs > 0 {
            (
                (stats.intra_pairs as f64 * 100.0) / total_pairs as f64,
                (stats.inter_pairs as f64 * 100.0) / total_pairs as f64,
            )
        } else {
            warn!("input BAM file '{}' contained 0 valid interaction pairs.", fname);
            (0.0, 0.0)
        };

        Self {
            fname,
            bin_size,
            filters,
            total_pairs,
            intra_pairs: stats.intra_pairs,
            inter_pairs: stats.inter_pairs,
            populated_bins: stats.bin_counts.len(),
            intra_pct,
            inter_pct,
        }
    }

    fn metrics(&self) -> Vec<(&'static str, String)> {
        vec![
            ("File name", self.fname.to_string()),
            ("Bin size (bp)", num_to_str(self.bin_size as u64)),
            ("Min MAPQ threshold", self.filters.min_mapq.to_string()),
            (
                "Max NM threshold",
                self.filters.max_nm.map_or_else(|| "None".to_string(), |v| v.to_string()),
            ),
            (
                "Min AS threshold",
                self.filters.min_as.map_or_else(|| "None".to_string(), |v| v.to_string()),
            ),
            ("Total read pairs", num_to_str(self.total_pairs as u64)),
            (
                "Intra-chromosomal pairs",
                format!("{} ({:.2}%)", num_to_str(self.intra_pairs as u64), self.intra_pct),
            ),
            (
                "Inter-chromosomal pairs",
                format!("{} ({:.2}%)", num_to_str(self.inter_pairs as u64), self.inter_pct),
            ),
            ("Populated matrix bins", num_to_str(self.populated_bins as u64)),
        ]
    }

    fn log_summary(&self) {
        info!("=================== MATRIX SUMMARY ===================");
        for (label, val) in self.metrics() {
            info!("{:<width$} : {}", label, val, width = DEFAULT_COLUMN_WIDTH);
        }
        info!("======================================================");
    }

    fn write_tsv(&self, bin_counts: &HashMap<ContactBinKey, u64>) -> Result<()> {
        let output_path = format!("{}.{}bp.matrix.tsv", self.fname, self.bin_size);
        let file = create_file(&output_path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "chrom1\tbin1\tchrom2\tbin2\tcount")?;

        let mut keys: Vec<&ContactBinKey> = bin_counts.keys().collect();
        keys.sort_by(|a, b| {
            a.chrom1
                .cmp(&b.chrom1)
                .then(a.bin1.cmp(&b.bin1))
                .then(a.chrom2.cmp(&b.chrom2))
                .then(a.bin2.cmp(&b.bin2))
        });

        for key in keys {
            let count = bin_counts.get(key).unwrap_or(&0);
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                key.chrom1, key.bin1, key.chrom2, key.bin2, count
            )?;
        }

        writer.flush()?;
        info!("saved sparse contact matrix to: {}", output_path);
        Ok(())
    }
}