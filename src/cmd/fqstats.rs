use crate::data::fastq::{FastqReader, FastqRecord};
use crate::error::{AppError, Result};
use crate::utils::helpers::create_file;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::io::{BufWriter, Write};

/// Intermediate stats collected per chunk/record across threads
#[derive(Default)]
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

    if nthreads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(nthreads)
            .build_global();
    }

    let reader = FastqReader::from_path(fname)?;

    const CHUNK_SIZE: usize = 10_000;
    let mut chunk_buffer: Vec<FastqRecord> = Vec::with_capacity(CHUNK_SIZE);
    let mut combined_stats = ChunkStats::new();

    for record_res in reader {
        let rec = record_res?;
        chunk_buffer.push(rec);

        if chunk_buffer.len() >= CHUNK_SIZE {
            let stats = process_chunk(&chunk_buffer);
            combined_stats = combined_stats.merge(stats);
            chunk_buffer.clear();
        }
    }

    if !chunk_buffer.is_empty() {
        let stats = process_chunk(&chunk_buffer);
        combined_stats = combined_stats.merge(stats);
    }

    let total_reads = combined_stats.total_reads;
    let total_bases = combined_stats.total_bases;

    let (min_len, max_len) = if total_reads == 0 {
        warn!("input file '{}' contained 0 records.", fname);
        (0, 0)
    } else {
        (combined_stats.min_length, combined_stats.max_length)
    };

    let offset: u8 = 33;

    let (min_qual, max_qual) = if total_reads == 0 {
        (offset, offset)
    } else {
        (combined_stats.min_qual_char, combined_stats.max_qual_char)
    };

    let min_phred = min_qual.saturating_sub(offset);
    let max_phred = max_qual.saturating_sub(offset);

    let pct_q20 = if total_bases > 0 {
        (combined_stats.bases_ge_q20 as f64 * 100.0) / total_bases as f64
    } else {
        0.0
    };

    let pct_q30 = if total_bases > 0 {
        (combined_stats.bases_ge_q30 as f64 * 100.0) / total_bases as f64
    } else {
        0.0
    };

    let avg_len = if total_reads > 0 {
        total_bases as f64 / total_reads as f64
    } else {
        0.0
    };

    let gc_count = combined_stats.count_g + combined_stats.count_c;
    let gc_pct = if total_bases > 0 {
        (gc_count as f64 * 100.0) / total_bases as f64
    } else {
        0.0
    };

    // Output strictly matching your requested format
    info!("File name : {}", fname);
    info!("Total bases : {}", total_bases);
    info!("Total reads : {}", total_reads);
    info!("% bases >=Q20 : {:.3}", pct_q20);
    info!("% bases >=Q30 : {:.3}", pct_q30);
    info!("Average read length : {:.3}", avg_len);
    info!("Read length range : {} .. {}", min_len, max_len);
    info!("Quality range : {} .. {}", min_qual, max_qual);
    info!("Phred range : {} .. {}", min_phred, max_phred);
    info!("Offset : {}", offset);
    info!("A : {}", combined_stats.count_a);
    info!("T : {}", combined_stats.count_t);
    info!("G : {}", combined_stats.count_g);
    info!("C : {}", combined_stats.count_c);
    info!("N : {}", combined_stats.count_n);
    info!("percent G-C content : {:.3}", gc_pct);

    let tsv_path = format!("{}.summary.tsv", fname);
    let file = create_file(&tsv_path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "File name\t{}", fname)?;
    writeln!(writer, "Total bases\t{}", total_bases)?;
    writeln!(writer, "Total reads\t{}", total_reads)?;
    writeln!(writer, "% bases >=Q20\t{:.3}", pct_q20)?;
    writeln!(writer, "% bases >=Q30\t{:.3}", pct_q30)?;
    writeln!(writer, "Average read length\t{:.3}", avg_len)?;
    writeln!(writer, "Read length range\t{} .. {}", min_len, max_len)?;
    writeln!(writer, "Quality range\t{} .. {}", min_qual, max_qual)?;
    writeln!(writer, "Phred range\t{} .. {}", min_phred, max_phred)?;
    writeln!(writer, "Offset\t{}", offset)?;
    writeln!(writer, "A\t{}", combined_stats.count_a)?;
    writeln!(writer, "T\t{}", combined_stats.count_t)?;
    writeln!(writer, "G\t{}", combined_stats.count_g)?;
    writeln!(writer, "C\t{}", combined_stats.count_c)?;
    writeln!(writer, "N\t{}", combined_stats.count_n)?;
    writeln!(writer, "percent G-C content\t{:.3}", gc_pct)?;

    writer.flush()?;
    debug!("saved summary metrics to TSV report: {}", tsv_path);

    Ok(())
}

fn process_chunk(records: &[FastqRecord]) -> ChunkStats {
    records
        .par_iter()
        .map(|rec| {
            let mut local = ChunkStats::new();
            let len = rec.len();

            local.total_reads = 1;
            local.total_bases = len;
            local.min_length = len;
            local.max_length = len;

            for &base in &rec.sequence {
                match base.to_ascii_uppercase() {
                    b'A' => local.count_a += 1,
                    b'T' => local.count_t += 1,
                    b'G' => local.count_g += 1,
                    b'C' => local.count_c += 1,
                    b'N' => local.count_n += 1,
                    _ => {}
                }
            }

            for &q in &rec.quality {
                local.min_qual_char = local.min_qual_char.min(q);
                local.max_qual_char = local.max_qual_char.max(q);

                let phred_score = q.saturating_sub(33);
                if phred_score >= 20 {
                    local.bases_ge_q20 += 1;
                }
                if phred_score >= 30 {
                    local.bases_ge_q30 += 1;
                }
            }

            local
        })
        .reduce(ChunkStats::new, |a, b| a.merge(b))
}