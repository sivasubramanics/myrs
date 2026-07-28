use crate::data::fasta::{FastaReader, FastaRecord};
use crate::data::kmc::Kmc;
use crate::data::kmer::Kmer;
use crate::error::{AppError, Result};
use crate::utils::helpers::create_file;
use crossbeam_channel::{bounded, Receiver, Sender};
use log::info;
use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
use fxhash::FxHashSet;

/// Result structure for a single processed FASTA sequence
struct RecordKmerStats {
    ref_name: String,
    ref_len: usize,
    ref_total_n: usize,
    ref_unique_n: usize,
    query_uniq_n: usize,
    query_total_n: u64,
    mean_ref: f64,
    mean_query: f64,
}

pub fn run(
    reference: &str,
    kmc_prefix: &str,
    output: &str,
    in_memory: bool,
    nthreads: usize,
) -> Result<()> {
    // 1. Clean KMC database path
    let db_base_path = kmc_prefix
        .strip_suffix(".kmc_pre")
        .or_else(|| kmc_prefix.strip_suffix(".kmc_suf"))
        .unwrap_or(kmc_prefix);

    let out_path = PathBuf::from(output);

    let kmc = Arc::new(Kmc::new(db_base_path, in_memory)?);
    let k = kmc.kmer_length();

    info!(
        "starting k-mer profiling on FASTA: '{}' (threads: {})",
        reference, nthreads
    );

    // Dedicate 1 thread for I/O reading, remainder for worker processing
    let worker_threads = if nthreads > 1 { nthreads - 1 } else { 1 };

    // Bounded channels keep RAM bounded when parsing huge FASTA files
    let (tx_records, rx_records): (Sender<FastaRecord>, Receiver<FastaRecord>) = bounded(100);
    let (tx_stats, rx_stats): (Sender<RecordKmerStats>, Receiver<RecordKmerStats>) = bounded(100);

    // -------------------------------------------------------------
    // 1. PRODUCER THREAD (FASTA I/O)
    // -------------------------------------------------------------
    let ref_path = reference.to_string();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let fasta_reader = FastaReader::from_path(&ref_path)?;
        for record_res in fasta_reader {
            let record = record_res?;
            if tx_records.send(record).is_err() {
                break; // Receiver disconnected
            }
        }
        Ok(())
    });

    // -------------------------------------------------------------
    // 2. WORKER THREADS (CPU / KMC Lookups)
    // -------------------------------------------------------------
    let mut worker_handles = Vec::with_capacity(worker_threads);

    for _ in 0..worker_threads {
        let rx = rx_records.clone();
        let tx = tx_stats.clone();
        let kmc_clone = Arc::clone(&kmc);

        let handle = thread::spawn(move || {
            while let Ok(record) = rx.recv() {
                let stats = process_fasta_record(&record, &kmc_clone, k);
                if tx.send(stats).is_err() {
                    break;
                }
            }
        });

        worker_handles.push(handle);
    }

    // Drop unused sender/receiver ends in the main thread so channels close naturally
    drop(rx_records);
    drop(tx_stats);

    // -------------------------------------------------------------
    // 3. AGGREGATOR & WRITER (Main Thread)
    // -------------------------------------------------------------
    info!("writing summary statistics to: '{}'", out_path.display());
    let file = create_file(&out_path)?;
    let mut writer = BufWriter::new(file);

    // Write TSV Header matching specified columns
    writeln!(
        writer,
        "ref\tref_len\tref_unique_n\tref_total_n\tquery_uniq_n\tquery_total_n\tmean_ref\tmean_query"
    )
        .map_err(AppError::Io)?;

    while let Ok(stats) = rx_stats.recv() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            stats.ref_name,
            stats.ref_len,
            stats.ref_unique_n,
            stats.ref_total_n,
            stats.query_uniq_n,
            stats.query_total_n,
            stats.mean_ref,
            stats.mean_query
        )
            .map_err(AppError::Io)?;
    }

    // Ensure reader finished cleanly without I/O errors
    reader_handle.join().unwrap()?;
    for handle in worker_handles {
        handle.join().unwrap();
    }

    writer.flush().map_err(AppError::Io)?;
    info!("successfully calculated k-mer statistics for '{}'.", reference);

    Ok(())
}


/// 2-bit encodes an ASCII DNA slice into a fixed [u64; 4] stack buffer (up to k=128).
/// Returns [u64; 4] representing the k-mer.
#[inline]
fn pack_kmer_2bit(kmer_bytes: &[u8]) -> Option<[u64; 4]> {
    let mut packed = [0u64; 4];
    for (i, &base) in kmer_bytes.iter().enumerate() {
        let val = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None, // Reject k-mers containing 'N' or non-standard bases
        };
        let word_idx = i / 32;
        let bit_shift = (i % 32) * 2;
        packed[word_idx] |= val << bit_shift;
    }
    Some(packed)
}

/// Compute KMC k-mer statistics for a single FASTA sequence record
#[inline]
fn process_fasta_record(record: &FastaRecord, kmc: &Kmc, k: usize) -> RecordKmerStats {
    if record.len() < k {
        return RecordKmerStats {
            ref_name: record.name.clone(),
            ref_len: record.len(),
            ref_unique_n: 0,
            ref_total_n: 0,
            query_uniq_n: 0,
            query_total_n: 0,
            mean_ref: 0.0,
            mean_query: 0.0,
        };
    }

    let mut total_kmers_ref = 0usize;

    // Fast non-cryptographic set storing zero-allocation [u64; 4] keys.
    // Perfectly supports any k up to 128!
    let mut unique_ref_kmers: FxHashSet<[u64; 4]> = FxHashSet::default();

    // Pre-allocate set capacity to avoid dynamic re-hashes during chromosome scanning
    unique_ref_kmers.reserve(record.len().saturating_sub(k - 1));

    let mut obs_unique_kmers = 0usize;
    let mut obs_kmer_total_count = 0u64;

    // Single-pass sliding window
    for kmer_ref in record.canonical_kmers(k) {
        total_kmers_ref += 1;

        info!("Processing k-mer: {:?}", std::str::from_utf8(kmer_ref.sequence()).unwrap_or("Invalid UTF-8"));

        // 1. Get raw ASCII bytes for canonical k-mer
        let bytes = kmer_ref.canonical_bytes();

        // 2. Pack into a zero-allocation 256-bit stack value
        let packed_key = pack_kmer_2bit(&bytes);

        // 3. Insert into FxHashSet. Returns true ONLY if newly encountered.
        if let Some(packed_key) = pack_kmer_2bit(&bytes) {
            if unique_ref_kmers.insert(packed_key) {
                let kmer = Kmer::new(&bytes);
                let count = kmc.get_kmer_count(&kmer) as u64;
                info!("Checking k-mer: {:?}, Count: {}", std::str::from_utf8(kmer.sequence()).unwrap_or("Invalid UTF-8"), count);
                if count > 0 {
                    info!("Observed k-mer: {:?}, Count: {}", std::str::from_utf8(kmer.sequence()).unwrap_or("Invalid UTF-8"), count);
                    obs_unique_kmers += 1;
                    obs_kmer_total_count += count;
                }
            }
        }
    }

    let n_unique_ref_kmers = unique_ref_kmers.len();

    let mean_ref = if n_unique_ref_kmers > 0 {
        obs_kmer_total_count as f64 / n_unique_ref_kmers as f64
    } else {
        0.0
    };

    let mean_query = if obs_unique_kmers > 0 {
        obs_kmer_total_count as f64 / obs_unique_kmers as f64
    } else {
        0.0
    };

    RecordKmerStats {
        ref_name: record.name.clone(),
        ref_len: record.len(),
        ref_unique_n: n_unique_ref_kmers,
        ref_total_n: total_kmers_ref,
        query_uniq_n: obs_unique_kmers,
        query_total_n: obs_kmer_total_count,
        mean_ref,
        mean_query,
    }
}