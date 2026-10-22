use crate::data::bam::{BamRecord, NativeBamReader, NativeBamWriter};
use crate::error::Result;
use crate::utils::helpers::create_file;
use log::{info, warn};
use std::io::{BufWriter, Write};

pub struct LinkOptions {
    pub filtered_bam: Option<String>,
    pub rejected_bam: Option<String>,
    pub only_primary: bool,
    pub max_nm: u8,
    pub min_mapq: u8,
}

pub fn run(
    input_path: &str,
    output_path: &str,
    _threads: usize,
    options: LinkOptions,
) -> Result<()> {
    let mut reader = NativeBamReader::from_path(input_path)?;
    let ref_names: Vec<String> = reader.header.references.iter().map(|r| r.name.clone()).collect();

    // Optional BAM writer for passing records
    let mut bam_writer = match options.filtered_bam {
        Some(ref path) => Some(NativeBamWriter::from_path(path, &reader.header)?),
        None => None,
    };

    // Optional BAM writer for rejected/filtered-out records
    let mut rejected_bam_writer = match options.rejected_bam {
        Some(ref path) => Some(NativeBamWriter::from_path(path, &reader.header)?),
        None => None,
    };

    let file = create_file(output_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "read_1_ref\tread_1_pos\tread_1_strand\tread_2_ref\tread_2_pos\tread_2_strand")?;

    let mut current = BamRecord::default();
    let mut pending: Option<BamRecord> = None;
    let mut records_seen: u64 = 0;
    let mut links_written: u64 = 0;
    let mut records_rejected: u64 = 0;

    while reader.read_next(&mut current)? {
        records_seen += 1;
        if options.only_primary && (current.is_secondary() || current.is_supplementary()) {
            if let Some(ref mut rbw) = rejected_bam_writer {
                rbw.write_record(&current)?;
            }
            records_rejected += 1;
            continue;
        }

        if let Some(prev) = pending.take() {
            if prev.name == current.name {
                let (r1, r2) = if prev.is_first_in_pair() && current.is_second_in_pair() {
                    (prev, current)
                } else if prev.is_second_in_pair() && current.is_first_in_pair() {
                    (current, prev)
                } else {
                    pending = Some(current);
                    current = BamRecord::default();
                    continue;
                };
                current = BamRecord::default();

                let i1 = r1.ref_id as usize;
                let i2 = r2.ref_id as usize;
                let valid_coords = i1 < ref_names.len() && i2 < ref_names.len() && r1.pos >= 0 && r2.pos >= 0;

                if !valid_pair(&r1, &r2, options.max_nm, options.min_mapq) || !valid_coords {
                    if let Some(ref mut rbw) = rejected_bam_writer {
                        rbw.write_record(&r1)?;
                        rbw.write_record(&r2)?;
                    }
                    records_rejected += 2;
                    continue;
                }

                if let Some(ref mut bw) = bam_writer {
                    bw.write_record(&r1)?;
                    bw.write_record(&r2)?;
                }

                let strand1 = if r1.flag & 0x10 != 0 { '-' } else { '+' };
                let strand2 = if r2.flag & 0x10 != 0 { '-' } else { '+' };
                writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}", ref_names[i1], r1.pos, strand1, ref_names[i2], r2.pos, strand2)?;
                links_written += 1;
            } else {
                pending = Some(current);
                current = BamRecord::default();
            }
        } else {
            pending = Some(current);
            current = BamRecord::default();
        }
    }

    writer.flush()?;
    if let Some(ref mut bw) = bam_writer {
        bw.write_eof_marker()?;
        bw.flush()?;
    }
    if let Some(ref mut rbw) = rejected_bam_writer {
        rbw.write_eof_marker()?;
        rbw.flush()?;
    }

    if links_written == 0 {
        warn!("no read links passed filtering; check BAM ordering and thresholds");
    }
    info!(
        "read {} BAM records; wrote {} links to {}; rejected {} records",
        records_seen, links_written, output_path, records_rejected
    );
    Ok(())
}

#[inline]
fn valid_pair(r1: &BamRecord, r2: &BamRecord, max_nm: u8, min_mapq: u8) -> bool {
    if !r1.is_paired() || !r2.is_paired()
        || r1.is_unmapped() || r2.is_unmapped()
        || r1.ref_id < 0 || r2.ref_id < 0
        || r1.mapq < min_mapq || r2.mapq < min_mapq
    {
        return false;
    }

    if r1.nm.map_or(true, |nm| nm > max_nm as i32)
        || r2.nm.map_or(true, |nm| nm > max_nm as i32)
    {
        return false;
    }
    true
}
