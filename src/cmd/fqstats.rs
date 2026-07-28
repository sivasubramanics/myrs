use crate::data::fastq::{FastqReader, FastqRecord};
use crate::error::{AppError, Result};
use crate::utils::helpers::{create_file, num_to_str};
use log::{debug, info, warn};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use crate::utils::DEFAULT_COLUMN_WIDTH;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::thread;


const PHRED_OFFSET: u8 = 33;
const CHUNK_SIZE: usize = 10_000;

// Pre-computed lookup table for fast nucleotide indexing
// Maps ASCII byte directly to 0=A, 1=T, 2=G, 3=C, 4=N, 255=Other
const BASE_LUT: [u8; 256] = {
    let mut lut = [255u8; 256];
    lut[b'A' as usize] = 0; lut[b'a' as usize] = 0;
    lut[b'T' as usize] = 1; lut[b't' as usize] = 1;
    lut[b'G' as usize] = 2; lut[b'g' as usize] = 2;
    lut[b'C' as usize] = 3; lut[b'c' as usize] = 3;
    lut[b'N' as usize] = 4; lut[b'n' as usize] = 4;
    lut
};

#[derive(Default, Debug, Clone)]
struct ChunkStats {
    total_reads: usize,
    total_bases: usize,
    min_length: usize,
    max_length: usize,
    bases_ge_q20: usize,
    bases_ge_q30: usize,
    min_qual_char: u8,
    max_qual_char: u8,
    count_a: usize,
    count_t: usize,
    count_g: usize,
    count_c: usize,
    count_n: usize,
}

impl ChunkStats {
    fn new() -> Self {
        Self {
            min_length: usize::MAX,
            min_qual_char: u8::MAX,
            ..Default::default()
        }
    }

    #[inline]
    fn merge(mut self, other: Self) -> Self {
        if other.total_reads == 0 {
            return self;
        }
        if self.total_reads == 0 {
            return other;
        }

        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;
        self.min_length = self.min_length.min(other.min_length);
        self.max_length = self.max_length.max(other.max_length);

        self.bases_ge_q20 += other.bases_ge_q20;
        self.bases_ge_q30 += other.bases_ge_q30;
        self.min_qual_char = self.min_qual_char.min(other.min_qual_char);
        self.max_qual_char = self.max_qual_char.max(other.max_qual_char);

        self.count_a += other.count_a;
        self.count_t += other.count_t;
        self.count_g += other.count_g;
        self.count_c += other.count_c;
        self.count_n += other.count_n;

        self
    }
}


pub fn run(fname: &str, nthreads: usize) -> Result<()> {
    info!("starting summary on file: {} (threads: {})", fname, nthreads);

    // Dedicate 1 thread for I/O reading, remainder for worker processing
    let worker_threads = if nthreads > 1 { nthreads - 1 } else { 1 };

    // Bounded channels prevent reading the entire file into RAM if workers lag
    let (tx_chunks, rx_chunks): (Sender<Vec<FastqRecord>>, Receiver<Vec<FastqRecord>>) = bounded(100);
    let (tx_stats, rx_stats): (Sender<ChunkStats>, Receiver<ChunkStats>) = bounded(100);

    // -------------------------------------------------------------
    // 1. PRODUCER THREAD (I/O)
    // -------------------------------------------------------------
    let fname_clone = fname.to_string();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let reader = FastqReader::from_path(&fname_clone)?;
        let mut chunk_buffer = Vec::with_capacity(CHUNK_SIZE);

        for record_res in reader {
            let rec = record_res?;
            chunk_buffer.push(rec);

            if chunk_buffer.len() >= CHUNK_SIZE {
                let full_chunk = std::mem::take(&mut chunk_buffer);
                if tx_chunks.send(full_chunk).is_err() {
                    break; // Receiver disconnected
                }
                chunk_buffer.reserve(CHUNK_SIZE);
            }
        }

        if !chunk_buffer.is_empty() {
            let _ = tx_chunks.send(chunk_buffer);
        }

        // `tx_chunks` is dropped automatically when this thread terminates!
        Ok(())
    });

    // -------------------------------------------------------------
    // 2. CONSUMER THREADS (CPU Computation)
    // -------------------------------------------------------------
    let mut worker_handles = Vec::with_capacity(worker_threads);

    for _ in 0..worker_threads {
        let rx = rx_chunks.clone();
        let tx = tx_stats.clone();

        let handle = thread::spawn(move || {
            let mut local_stats = ChunkStats::new();

            while let Ok(chunk) = rx.recv() {
                for rec in chunk {
                    process_single_record(&rec, &mut local_stats);
                }
            }

            // Send worker results back to main thread
            let _ = tx.send(local_stats);
        });

        worker_handles.push(handle);
    }

    drop(rx_chunks);
    drop(tx_stats);

    // -------------------------------------------------------------
    // 3. AGGREGATOR (Main Thread)
    // -------------------------------------------------------------
    let mut combined_stats = ChunkStats::new();
    while let Ok(stats) = rx_stats.recv() {
        combined_stats = combined_stats.merge(stats);
    }

    // Ensure reader finished cleanly without I/O errors
    reader_handle.join().unwrap()?;
    for handle in worker_handles {
        handle.join().unwrap();
    }

    let summary = CalculatedSummary::from_stats(fname, &combined_stats);
    summary.log_summary();
    summary.write_tsv()?;

    Ok(())
}

#[inline]
fn process_single_record(rec: &FastqRecord, stats: &mut ChunkStats) {
    let len = rec.sequence.len();
    stats.total_reads += 1;
    stats.total_bases += len;
    stats.min_length = stats.min_length.min(len);
    stats.max_length = stats.max_length.max(len);

    let q20_threshold = PHRED_OFFSET + 20;
    let q30_threshold = PHRED_OFFSET + 30;

    // Single pass over sequence
    for &b in &rec.sequence {
        match BASE_LUT[b as usize] {
            0 => stats.count_a += 1,
            1 => stats.count_t += 1,
            2 => stats.count_g += 1,
            3 => stats.count_c += 1,
            4 => stats.count_n += 1,
            _ => {}
        }
    }

    // Single pass over quality
    for &q in &rec.quality {
        stats.min_qual_char = stats.min_qual_char.min(q);
        stats.max_qual_char = stats.max_qual_char.max(q);
        if q >= q20_threshold { stats.bases_ge_q20 += 1; }
        if q >= q30_threshold { stats.bases_ge_q30 += 1; }
    }
}

/// Optimized multi-threaded chunk processor using lookup tables
fn process_chunk(records: &[FastqRecord]) -> ChunkStats {
    records
        .par_iter()
        .map(|rec| {
            let mut local = ChunkStats::new();
            let len = rec.sequence.len();

            local.total_reads = 1;
            local.total_bases = len;
            local.min_length = len;
            local.max_length = len;

            let mut counts = [0usize; 5];
            let q20_threshold = PHRED_OFFSET + 20;
            let q30_threshold = PHRED_OFFSET + 30;

            // Tight loop for sequence base counting via LUT
            for &b in &rec.sequence {
                let idx = BASE_LUT[b as usize];
                if idx < 5 {
                    counts[idx as usize] += 1;
                }
            }

            local.count_a = counts[0];
            local.count_t = counts[1];
            local.count_g = counts[2];
            local.count_c = counts[3];
            local.count_n = counts[4];

            // Tight loop for quality scores
            for &q in &rec.quality {
                local.min_qual_char = local.min_qual_char.min(q);
                local.max_qual_char = local.max_qual_char.max(q);

                if q >= q20_threshold {
                    local.bases_ge_q20 += 1;
                }
                if q >= q30_threshold {
                    local.bases_ge_q30 += 1;
                }
            }

            local
        })
        .reduce(ChunkStats::new, |a, b| a.merge(b))
}

struct CalculatedSummary<'a> {
    fname: &'a str,
    total_bases: usize,
    total_reads: usize,
    pct_q20: f64,
    pct_q30: f64,
    avg_len: f64,
    min_len: usize,
    max_len: usize,
    min_qual: u8,
    max_qual: u8,
    min_phred: u8,
    max_phred: u8,
    count_a: usize,
    count_t: usize,
    count_g: usize,
    count_c: usize,
    count_n: usize,
    gc_pct: f64,
}

impl<'a> CalculatedSummary<'a> {
    fn from_stats(fname: &'a str, stats: &ChunkStats) -> Self {
        let total_reads = stats.total_reads;
        let total_bases = stats.total_bases;

        let (min_len, max_len) = if total_reads == 0 {
            warn!("input file '{}' contained 0 records.", fname);
            (0, 0)
        } else {
            (stats.min_length, stats.max_length)
        };

        let (min_qual, max_qual) = if total_reads == 0 {
            (PHRED_OFFSET, PHRED_OFFSET)
        } else {
            (stats.min_qual_char, stats.max_qual_char)
        };

        let min_phred = min_qual.saturating_sub(PHRED_OFFSET);
        let max_phred = max_qual.saturating_sub(PHRED_OFFSET);

        let pct_q20 = if total_bases > 0 {
            (stats.bases_ge_q20 as f64 * 100.0) / total_bases as f64
        } else {
            0.0
        };

        let pct_q30 = if total_bases > 0 {
            (stats.bases_ge_q30 as f64 * 100.0) / total_bases as f64
        } else {
            0.0
        };

        let avg_len = if total_reads > 0 {
            total_bases as f64 / total_reads as f64
        } else {
            0.0
        };

        let gc_count = stats.count_g + stats.count_c;
        let gc_pct = if total_bases > 0 {
            (gc_count as f64 * 100.0) / total_bases as f64
        } else {
            0.0
        };

        Self {
            fname,
            total_bases,
            total_reads,
            pct_q20,
            pct_q30,
            avg_len,
            min_len,
            max_len,
            min_qual,
            max_qual,
            min_phred,
            max_phred,
            count_a: stats.count_a,
            count_t: stats.count_t,
            count_g: stats.count_g,
            count_c: stats.count_c,
            count_n: stats.count_n,
            gc_pct,
        }
    }

    /// Zero-allocation key-value pair iterator
    fn metrics(&self) -> Vec<(&'static str, String)> {
        vec![
            ("File name", self.fname.to_string()),
            ("Total bases", num_to_str(self.total_bases as u64)),
            ("Total reads", num_to_str(self.total_reads as u64)),
            ("% bases >=Q20", format!("{:.3}", self.pct_q20)),
            ("% bases >=Q30", format!("{:.3}", self.pct_q30)),
            ("Average read length", format!("{:.3}", self.avg_len)),
            ("Read length range", format!("{} .. {}", self.min_len, self.max_len)),
            ("Quality range", format!("{} .. {}", self.min_qual, self.max_qual)),
            ("Phred range", format!("{} .. {}", self.min_phred, self.max_phred)),
            ("Offset", PHRED_OFFSET.to_string()),
            ("A", num_to_str(self.count_a as u64)),
            ("T", num_to_str(self.count_t as u64)),
            ("G", num_to_str(self.count_g as u64)),
            ("C", num_to_str(self.count_c as u64)),
            ("N", num_to_str(self.count_n as u64)),
            ("percent G-C content", format!("{:.3}", self.gc_pct)),
        ]
    }

    fn log_summary(&self) {
        info!("====================== SUMMARY ======================");
        for (label, val) in self.metrics() {
            info!("{:<width$} : {}", label, val, width = DEFAULT_COLUMN_WIDTH);
        }
        info!("=====================================================");
    }

    fn write_tsv(&self) -> Result<()> {
        let tsv_path = format!("{}.summary.tsv", self.fname);
        let file = create_file(&tsv_path)?;
        let mut writer = BufWriter::new(file);

        for (label, val) in self.metrics() {
            writeln!(writer, "{}\t{}", label, val)?;
        }

        writer.flush()?;
        debug!("saved summary metrics to TSV report: {}", tsv_path);
        Ok(())
    }
}