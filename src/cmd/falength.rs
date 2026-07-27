use std::io::{BufWriter, Write};
use log::info;
use crate::data::fasta::{FastaReader, FastaRecord};
use crate::utils::helpers::create_file;

pub fn run(fname: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = FastaReader::from_path(fname)?;

    let outfile = create_file(&format!("{}.len", fname))?;
    let mut writer = BufWriter::new(outfile);

    // Reuse a single record object across every iteration
    let mut record = FastaRecord::default();

    while reader.read_next(&mut record)? {
        writeln!(writer, "{}\t{}", record.name, record.len())?;
    }

    // Flush the writer to ensure all data is saved
    writer.flush()?;
    info!("length information written to '{}.len'", fname);

    Ok(())
}