/// Default buffer size for buffered I/O readers/writers (512 KB).
pub const DEFAULT_BUF_SIZE: usize = 512 * 1024;

/// Default line wrap width for FASTA output sequence lines.
pub const DEFAULT_FOLD_WIDTH: usize = 80;

/// Standard FASTA file extensions recognized across bioinformatics tools.
pub const VALID_FASTA_EXTENSIONS: &[&str] = &["fasta", "fa", "fna", "faa", "ffn", "frn"];

/// Standard FASTQ file extensions.
pub const VALID_FASTQ_EXTENSIONS: &[&str] = &["fastq", "fq"];

/// column width
pub const DEFAULT_COLUMN_WIDTH: usize = 30;

// 256-byte lookup table mapping every ASCII character to its reverse complement.
/// Preserves upper/lowercase and supports standard IUPAC ambiguity codes.
pub const COMPLEMENT: [u8; 256] = {
    let mut table = [b'N'; 256];

    // Uppercase standard
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'G' as usize] = b'C';
    table[b'C' as usize] = b'G';
    table[b'U' as usize] = b'A';

    // Lowercase standard
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'g' as usize] = b'c';
    table[b'c' as usize] = b'g';
    table[b'u' as usize] = b'a';

    // IUPAC Ambiguity Codes (Uppercase)
    table[b'R' as usize] = b'Y'; table[b'Y' as usize] = b'R';
    table[b'S' as usize] = b'S'; table[b'W' as usize] = b'W';
    table[b'K' as usize] = b'M'; table[b'M' as usize] = b'K';
    table[b'B' as usize] = b'V'; table[b'V' as usize] = b'B';
    table[b'D' as usize] = b'H'; table[b'H' as usize] = b'D';
    table[b'N' as usize] = b'N';

    // IUPAC Ambiguity Codes (Lowercase)
    table[b'r' as usize] = b'y'; table[b'y' as usize] = b'r';
    table[b's' as usize] = b's'; table[b'w' as usize] = b'w';
    table[b'k' as usize] = b'm'; table[b'm' as usize] = b'k';
    table[b'b' as usize] = b'v'; table[b'v' as usize] = b'b';
    table[b'd' as usize] = b'h'; table[b'h' as usize] = b'd';
    table[b'n' as usize] = b'n';

    table
};