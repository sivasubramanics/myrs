use crate::data::fasta::{FastaReader, FastaRecord};
use crate::data::kmc::Kmc;
use crate::error::{AppError, Result};
use crate::utils::helpers::{create_file, num_to_str};
use crossbeam_channel::{bounded, Receiver, Sender};
use fxhash::FxHashSet;
use log::info;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;


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

/// Convert a k-mer byte slice (ASCII A, C, G, T) into a bit-packed 2-bit u64
#[inline(always)]
fn encode_kmer(slice: &[u8]) -> Option<u64> {
    let mut val: u64 = 0;
    for &b in slice {
        let code = match b {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => return None, // Non-ACGT characters (e.g., 'N')
        };
        val = (val << 2) | code;
    }
    Some(val)
}

fn process_fasta_record(record: &FastaRecord, kmc: &Kmc, k: usize) -> RecordKmerStats {
    let ref_name = record.name.clone();
    let ref_len = record.len();

    // Zero-allocation hash set using packed u64 k-mers (for K <= 32)
    // Pre-allocate capacity to reduce re-hash allocations
    let estimated_kmers = ref_len.saturating_sub(k - 1);
    let mut seen_kmers: FxHashSet<u64> = FxHashSet::with_capacity_and_hasher(
        estimated_kmers,
        Default::default()
    );

    let mut query_total_n: u64 = 0;
    let mut query_uniq_n: usize = 0;
    let mut ref_total_n: usize = 0;

    for kmer_bytes in record.kmers(k, true) {
        ref_total_n += 1;

        // Query KMC count using byte slice
        let count = kmc.get_count(&kmer_bytes);
        if count > 0 {
            query_total_n += count as u64;

            if let Some(packed) = encode_kmer(&kmer_bytes) {
                if seen_kmers.insert(packed) {
                    query_uniq_n += 1;
                }
            }
        } else if let Some(packed) = encode_kmer(&kmer_bytes) {
            seen_kmers.insert(packed);
        }
    }

    let ref_unique_n = seen_kmers.len();

    let mean_ref = if ref_unique_n > 0 {
        query_total_n as f64 / ref_unique_n as f64
    } else {
        0.0
    };

    let mean_query = if ref_total_n > 0 {
        query_total_n as f64 / ref_total_n as f64
    } else {
        0.0
    };

    RecordKmerStats {
        ref_name,
        ref_len,
        ref_total_n,
        ref_unique_n,
        query_uniq_n,
        query_total_n,
        mean_ref,
        mean_query,
    }
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
    let (tx_records, rx_records): (Sender<FastaRecord>, Receiver<FastaRecord>) = bounded(1000);
    let (tx_stats, rx_stats): (Sender<RecordKmerStats>, Receiver<RecordKmerStats>) = bounded(1000);

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

    let mut processed_records: usize = 0;
    let mut stderr = io::stderr();

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

        processed_records += 1;

        // Print in-place progress update on stderr
        // Update frequency can be throttled (e.g., `if processed_records % 100 == 0`) for high throughput
        if processed_records % 1000 == 0 {
            eprint!("\rProcessed FASTA records: {}", num_to_str(processed_records as u64));
            let _ = stderr.flush();
        }
        let _ = stderr.flush();
    }

    // Clear the progress line after completion
    eprintln!("\rProcessed FASTA records: {} (done)", num_to_str(processed_records as u64));

    // Ensure reader finished cleanly without I/O errors
    reader_handle.join().unwrap()?;
    for handle in worker_handles {
        handle.join().unwrap();
    }

    writer.flush().map_err(AppError::Io)?;
    info!("successfully calculated k-mer statistics for '{}'.", reference);

    Ok(())
}

// use crate::data::fasta::{FastaReader, FastaRecord};
// use crate::data::kmc::Kmc;
// use crate::error::{AppError, Result};
// use crate::utils::helpers::{create_file, num_to_str};
// use fxhash::FxHashSet;
// use log::info;
// use rayon::prelude::*;
// use std::io::{self, BufWriter, Write};
// use std::path::PathBuf;
// use std::sync::atomic::{AtomicUsize, Ordering};
// use std::sync::{Arc, Mutex};
// use rayon::iter::ParallelBridge;
// use rayon::prelude::*;
//
// /// Result structure for a single processed FASTA sequence
// struct RecordKmerStats {
//     ref_name: String,
//     ref_len: usize,
//     ref_total_n: usize,
//     ref_unique_n: usize,
//     query_uniq_n: usize,
//     query_total_n: u64,
//     mean_ref: f64,
//     mean_query: f64,
// }
//
// /// Convert a k-mer byte slice (ASCII A, C, G, T) into a bit-packed 2-bit u64
// #[inline(always)]
// fn encode_kmer(slice: &[u8]) -> Option<u64> {
//     let mut val: u64 = 0;
//     for &b in slice {
//         let code = match b {
//             b'A' | b'a' => 0b00,
//             b'C' | b'c' => 0b01,
//             b'G' | b'g' => 0b10,
//             b'T' | b't' => 0b11,
//             _ => return None, // Non-ACGT characters (e.g., 'N')
//         };
//         val = (val << 2) | code;
//     }
//     Some(val)
// }
//
// fn process_fasta_record(record: &FastaRecord, kmc: &Kmc, k: usize) -> RecordKmerStats {
//     let ref_name = record.name.clone();
//     let ref_len = record.len();
//
//     // Zero-allocation hash set using packed u64 k-mers (for K <= 32)
//     let estimated_kmers = ref_len.saturating_sub(k - 1);
//     let mut seen_kmers: FxHashSet<u64> = FxHashSet::with_capacity_and_hasher(
//         estimated_kmers,
//         Default::default(),
//     );
//
//     let mut query_total_n: u64 = 0;
//     let mut query_uniq_n: usize = 0;
//     let mut ref_total_n: usize = 0;
//
//     for kmer_bytes in record.kmers(k, true) {
//         ref_total_n += 1;
//
//         // Query KMC count using byte slice
//         let count = kmc.get_count(&kmer_bytes);
//
//         if count > 0 {
//             query_total_n += count as u64;
//
//             if let Some(packed) = encode_kmer(&kmer_bytes) {
//                 if seen_kmers.insert(packed) {
//                     query_uniq_n += 1;
//                 }
//             }
//         } else if let Some(packed) = encode_kmer(&kmer_bytes) {
//             seen_kmers.insert(packed);
//         }
//     }
//
//     let ref_unique_n = seen_kmers.len();
//
//     let mean_ref = if ref_unique_n > 0 {
//         query_total_n as f64 / ref_unique_n as f64
//     } else {
//         0.0
//     };
//
//     let mean_query = if ref_total_n > 0 {
//         query_total_n as f64 / ref_total_n as f64
//     } else {
//         0.0
//     };
//
//     RecordKmerStats {
//         ref_name,
//         ref_len,
//         ref_total_n,
//         ref_unique_n,
//         query_uniq_n,
//         query_total_n,
//         mean_ref,
//         mean_query,
//     }
// }
//
// pub fn run(
//     reference: &str,
//     kmc_prefix: &str,
//     output: &str,
//     in_memory: bool,
//     nthreads: usize,
// ) -> Result<()> {
//     // 1. Clean KMC database path
//     let db_base_path = kmc_prefix
//         .strip_suffix(".kmc_pre")
//         .or_else(|| kmc_prefix.strip_suffix(".kmc_suf"))
//         .unwrap_or(kmc_prefix);
//
//     let out_path = PathBuf::from(output);
//
//     let kmc = Arc::new(Kmc::new(db_base_path, in_memory)?);
//     let k = kmc.kmer_length();
//
//     info!(
//         "starting k-mer profiling on FASTA: '{}' (threads: {})",
//         reference, nthreads
//     );
//
//     // Build dedicated Rayon thread pool matching nthreads parameter
//     let pool = rayon::ThreadPoolBuilder::new()
//         .num_threads(nthreads)
//         .build()
//         .map_err(|e| AppError::ThreadingError(e.to_string()))?;
//
//     // Output writer wrapped in Mutex for thread safety across Rayon workers
//     let file = create_file(&out_path)?;
//     let writer = Arc::new(Mutex::new(BufWriter::new(file)));
//
//     // Write TSV Header
//     {
//         let mut w = writer.lock().unwrap();
//         writeln!(
//             w,
//             "ref\tref_len\tref_unique_n\tref_total_n\tquery_uniq_n\tquery_total_n\tmean_ref\tmean_query"
//         )
//             .map_err(AppError::Io)?;
//     }
//
//     let mut fasta_reader = FastaReader::from_path(reference)?;
//     let processed_records = Arc::new(AtomicUsize::new(0));
//
//     // Execute parallel processing inside Rayon pool
//     pool.install(|| -> Result<()> {
//         // std::iter::from_fn creates an Iterator that implements Send cleanly
//         std::iter::from_fn(move || fasta_reader.next())
//             .par_bridge() // Parallelize single-threaded reader stream using work-stealing
//             .try_for_each(|record_res| -> Result<()> {
//                 let record = record_res?;
//                 let stats = process_fasta_record(&record, &kmc, k);
//
//                 // Write results incrementally
//                 {
//                     let mut w = writer.lock().unwrap();
//                     writeln!(
//                         w,
//                         "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
//                         stats.ref_name,
//                         stats.ref_len,
//                         stats.ref_unique_n,
//                         stats.ref_total_n,
//                         stats.query_uniq_n,
//                         stats.query_total_n,
//                         stats.mean_ref,
//                         stats.mean_query
//                     )
//                         .map_err(AppError::Io)?;
//                 }
//
//                 // Update progress atomically
//                 let count = processed_records.fetch_add(1, Ordering::Relaxed) + 1;
//                 if count % 1000 == 0 {
//                     let mut stderr = io::stderr();
//                     eprint!("\rProcessed FASTA records: {}", num_to_str(count as u64));
//                     let _ = stderr.flush();
//                 }
//
//                 Ok(())
//             })
//     })?;
//
//     let total = processed_records.load(Ordering::Relaxed);
//     eprintln!(
//         "\rProcessed FASTA records: {} (done)",
//         num_to_str(total as u64)
//     );
//
//     // Flush final output
//     let mut w = writer.lock().unwrap();
//     w.flush().map_err(AppError::Io)?;
//
//     info!(
//         "successfully calculated k-mer statistics for '{}'.",
//         reference
//     );
//
//     Ok(())
// }