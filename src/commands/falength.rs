use std::fs::File;
use std::io::{BufWriter, Write};

use crate::data::fasta::{FastaReader, FastaRecord};

pub fn run(fname: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = FastaReader::from_path(fname)?;

    let outfile = File::create(format!("{fname}.len"))?;
    let mut writer = BufWriter::new(outfile);

    // Reuse a single record object across every iteration
    let mut record = FastaRecord::default();

    while reader.read_next(&mut record)? {
        writeln!(writer, "{}\t{}", record.name, record.len())?;
    }

    // Flush the writer to ensure all data is saved
    writer.flush()?;

    Ok(())
}