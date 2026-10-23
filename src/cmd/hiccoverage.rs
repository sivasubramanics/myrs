use crate::error::Result;
use crate::utils::helpers::create_file;
use crossbeam_channel::{bounded, Receiver, Sender};
use log::info;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::thread;

const CHUNK_SIZE: usize = 50_000;

#[derive(Debug, Default, Clone)]
struct CoverageStats {
    total_lines: usize,
    valid_pairs: usize,
}

// Data structure to hold both local counters from worker threads
struct WorkerResult {
    coverage_counter: HashMap<(String, usize), usize>,
    matrix_counter: HashMap<(String, usize, String, usize), usize>,
    stats: CoverageStats,
}

pub fn run(
    links_path: &str,
    output_path: &str,
    matrix_output_path: Option<&str>,
    binsize: usize,
    nthreads: usize,
) -> Result<()> {
    if binsize == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "binsize must be greater than zero",
        )
            .into());
    }

    info!(
        "starting Hi-C processing from '{}' (binsize: {}, threads: {})",
        links_path, binsize, nthreads
    );

    let worker_threads = if nthreads > 1 { nthreads - 1 } else { 1 };

    let (tx_chunks, rx_chunks): (Sender<Vec<String>>, Receiver<Vec<String>>) = bounded(100);
    let (tx_results, rx_results): (
        Sender<std::result::Result<WorkerResult, String>>,
        Receiver<std::result::Result<WorkerResult, String>>,
    ) = bounded(100);

    // -------------------------------------------------------------
    // 1. PRODUCER THREAD (File I/O Streaming)
    // -------------------------------------------------------------
    let links_path_clone = links_path.to_string();
    let reader_handle = thread::spawn(move || -> io::Result<()> {
        let file = File::open(&links_path_clone)?;
        let reader = BufReader::with_capacity(64 * 1024, file);
        let mut chunk_buffer = Vec::with_capacity(CHUNK_SIZE);

        for line_res in reader.lines() {
            let line = line_res?;
            chunk_buffer.push(line);
            if chunk_buffer.len() >= CHUNK_SIZE {
                let full_chunk = std::mem::take(&mut chunk_buffer);
                if tx_chunks.send(full_chunk).is_err() {
                    break;
                }
                chunk_buffer.reserve(CHUNK_SIZE);
            }
        }

        if !chunk_buffer.is_empty() {
            let _ = tx_chunks.send(chunk_buffer);
        }

        Ok(())
    });

    // -------------------------------------------------------------
    // 2. WORKER THREADS (CPU Parallel Parsing & Binning)
    // -------------------------------------------------------------
    let matrix_enabled = matrix_output_path.is_some();
    let mut worker_handles = Vec::with_capacity(worker_threads);

    for _ in 0..worker_threads {
        let rx = rx_chunks.clone();
        let tx = tx_results.clone();

        let handle = thread::spawn(move || {
            let mut coverage_counter: HashMap<(String, usize), usize> = HashMap::default();
            let mut matrix_counter: HashMap<(String, usize, String, usize), usize> = HashMap::default();
            let mut stats = CoverageStats::default();

            while let Ok(chunk) = rx.recv() {
                for line in chunk {
                    stats.total_lines += 1;

                    // Field validation: expect exactly 6 fields
                    let mut parts = line.split('\t');
                    let r1_ref = parts.next();
                    let p1_str = parts.next();
                    let _r1_strand = parts.next();
                    let r2_ref = parts.next();
                    let p2_str = parts.next();
                    let r2_strand = parts.next();

                    if r1_ref.is_none()
                        || p1_str.is_none()
                        || _r1_strand.is_none()
                        || r2_ref.is_none()
                        || p2_str.is_none()
                        || r2_strand.is_none()
                        || parts.next().is_some()
                    {
                        let _ = tx.send(Err(format!(
                            "Invalid line format: expected exactly 6 fields. Line: {}",
                            line
                        )));
                        return;
                    }

                    let r1_ref = r1_ref.unwrap();
                    let p1_str = p1_str.unwrap();
                    let r2_ref = r2_ref.unwrap();
                    let p2_str = p2_str.unwrap();

                    let p1 = match p1_str.parse::<usize>() {
                        Ok(val) => val,
                        Err(_) => continue,
                    };
                    let p2 = match p2_str.parse::<usize>() {
                        Ok(val) => val,
                        Err(_) => continue,
                    };

                    // Compute bin indices for endpoints
                    let bin1 = (p1 / binsize) * binsize;
                    let bin2 = (p2 / binsize) * binsize;

                    // 1. Accumulate Contact Matrix if requested
                    if matrix_enabled {
                        // Canonicalize pair order for contact matrix if needed, or keep as-is
                        let (chrom_a, b_a, chrom_b, b_b) = if r1_ref < r2_ref || (r1_ref == r2_ref && bin1 <= bin2) {
                            (r1_ref.to_string(), bin1, r2_ref.to_string(), bin2)
                        } else {
                            (r2_ref.to_string(), bin2, r1_ref.to_string(), bin1)
                        };
                        *matrix_counter.entry((chrom_a, b_a, chrom_b, b_b)).or_insert(0) += 1;
                    }

                    // 2. Physical Coverage Filtering & Calculation (Intra-scaffold only)
                    if r1_ref != r2_ref {
                        continue;
                    }

                    let (p_start, p_end) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
                    let i = p_start / binsize;
                    let j = p_end / binsize;

                    if j - i <= 2 {
                        continue;
                    }

                    stats.valid_pairs += 1;

                    for bin_idx in (i + 1)..j {
                        *coverage_counter
                            .entry((r1_ref.to_string(), bin_idx * binsize))
                            .or_insert(0) += 1;
                    }
                }
            }

            let _ = tx.send(Ok(WorkerResult {
                coverage_counter,
                matrix_counter,
                stats,
            }));
        });

        worker_handles.push(handle);
    }

    drop(rx_chunks);
    drop(tx_results);

    // -------------------------------------------------------------
    // 3. AGGREGATOR & THREAD SHUTDOWN
    // -------------------------------------------------------------
    let mut global_coverage: HashMap<(String, usize), usize> = HashMap::default();
    let mut global_matrix: HashMap<(String, usize, String, usize), usize> = HashMap::default();
    let mut total_lines = 0;
    let mut valid_pairs = 0;
    let mut first_error: Option<String> = None;

    while let Ok(worker_res) = rx_results.recv() {
        match worker_res {
            Ok(res) => {
                total_lines += res.stats.total_lines;
                valid_pairs += res.stats.valid_pairs;

                for (key, count) in res.coverage_counter {
                    *global_coverage.entry(key).or_insert(0) += count;
                }
                for (key, count) in res.matrix_counter {
                    *global_matrix.entry(key).or_insert(0) += count;
                }
            }
            Err(err_msg) => {
                if first_error.is_none() {
                    first_error = Some(err_msg);
                }
            }
        }
    }

    let reader_res = reader_handle.join();
    match reader_res {
        Ok(Ok(())) => {}
        Ok(Err(e)) => {
            if first_error.is_none() {
                first_error = Some(e.to_string());
            }
        }
        Err(_) => {
            if first_error.is_none() {
                first_error = Some("Reader thread panicked".to_string());
            }
        }
    }

    for handle in worker_handles {
        if let Err(_) = handle.join() {
            if first_error.is_none() {
                first_error = Some("Worker thread panicked".to_string());
            }
        }
    }

    if let Some(err_msg) = first_error {
        return Err(io::Error::new(io::ErrorKind::InvalidData, err_msg).into());
    }

    info!(
        "Processed {} lines, found {} valid intra-scaffold pairs.",
        crate::utils::helpers::num_to_str(total_lines as u64),
        crate::utils::helpers::num_to_str(valid_pairs as u64)
    );

    // Write Physical Coverage Output
    info!("Writing coverage results to '{}'...", output_path);
    let mut sorted_coverage_keys: Vec<(String, usize)> = global_coverage.keys().cloned().collect();
    sorted_coverage_keys.sort_unstable();

    let out_file = create_file(output_path)?;
    let mut writer = BufWriter::new(out_file);
    for key in sorted_coverage_keys {
        let count = global_coverage[&key];
        writeln!(writer, "{}\t{}\t{}", key.0, key.1, count)?;
    }
    writer.flush()?;

    // Write Contact Matrix Output (if requested)
    if let Some(matrix_path) = matrix_output_path {
        info!("Writing contact matrix results to '{}'...", matrix_path);
        let mut sorted_matrix_keys: Vec<(String, usize, String, usize)> =
            global_matrix.keys().cloned().collect();
        sorted_matrix_keys.sort_unstable();

        let matrix_file = create_file(matrix_path)?;
        let mut matrix_writer = BufWriter::new(matrix_file);
        for key in sorted_matrix_keys {
            let count = global_matrix[&key];
            writeln!(
                matrix_writer,
                "{}\t{}\t{}\t{}\t{}",
                key.0, key.1, key.2, key.3, count
            )?;
        }
        matrix_writer.flush()?;
    }

    Ok(())
}