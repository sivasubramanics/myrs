# `myrs`

`myrs` is a simple project built to learn Rust. Many other tools already provide the same functionality.  

---

## Installation & Setup

To build `myrs` from source:

```bash
# Clone the repository
git clone https://github.com/sivasubramanics/myrs.git
cd myrs

# Build optimized release binary
cargo build --release

# Check the executable works
target/release/myrs --help

# Optional: Add to PATH or copy to group binaries directory
cp target/release/myrs ~/.local/bin/

```

---

## Core Commands & Usage

### 1. FASTA Utilities

#### `fa-stats`

Calculate summary statistics ($N_{50}$, total length, sequence count, GC%).

```bash
myrs fa-stats -f input.fasta.gz

```

* `-f, --fname <FILE>`: Input FASTA file path (supports `.gz`).



---

#### `fa-length`

Report sequence lengths for all records in a FASTA file.

```bash
myrs fa-length -f input.fasta 

```

* `-f, --fname <FILE>`: Input FASTA file path (supports `.gz`).



---

#### `fa-filter`

Filter FASTA records by sequence length and GC percentage.

```bash
myrs fa-filter -f input.fasta -m 1000 -M 50000 -g 40.0 > filtered.fasta

```

* `-f, --fname <FILE>`: Input FASTA file path (supports `.gz`).


* `-m, --min-len <INT>`: Minimum sequence length threshold (bp).


* `-M, --max-len <INT>`: Maximum sequence length threshold (bp).


* `-g, --min-gc <FLOAT>`: Minimum GC percentage threshold (e.g., `40.0` for 40%).



---

#### `fa-rename`

Rename FASTA headers using a 2-column TSV mapping file (`old_name\tnew_name`).

```bash
myrs fa-rename -f input.fasta -m mapping.tsv -o renamed.fasta --keep-old-id

```

* `-f, --fname <FILE>`: Input FASTA file path.


* `-m, --map <FILE>`: 2-column TSV mapping file (`old_name\tnew_name`).


* `-o, --out <FILE>`: Output FASTA file path (default: `<input>.renamed.fa`).


* `--keep-old-id`: Keep original sequence ID in the header description as `old_id="<name>"`.


* `-w, --width <INT>`: Line folding/wrapping width for output sequence lines.



---

### 2. Sequence Extraction & Indexing

#### `fa-fai`

Generate a `.faidx` index file for fast random access.

```bash
myrs fa-fai -f input.fasta

```

* `-f, --fname <FILE>`: Input FASTA file path (uncompressed).



---

#### `fa-one-record`

Rapidly extract a single sequence record from an indexed FASTA file.

```bash
myrs fa-one-record -f input.fasta -n "chr1" -o chr1.fasta

```

* `-f, --fname <FILE>`: Input FASTA file path.


* `-n, --name <STRING>`: Sequence ID or header name to extract.


* `-o, --out <FILE>`: Output file path (default: `<input>.<name>.fa`).


* `-w, --width <INT>`: Line folding/wrapping width for output sequence lines.



---

#### `fa-some-records`

Extract multiple sequence records by an ID list or list file.

```bash
# Using a file with one ID per line
myrs fa-some-records -f input.fasta -q headers.txt -o subset.fasta

# Using inline names (comma or space separated)
myrs fa-some-records -f input.fasta -n chr1,chr2,chr3 -o subset.fasta

```

* `-f, --fname <FILE>`: Input FASTA file path.


* `-q, --names-file <FILE>`: Text file containing sequence IDs to extract (one per line).


* `-n, --names <STRINGS>`: Sequence IDs to extract (comma or space separated).


* `-o, --out <FILE>`: Output file path (default: `<input>.some.fa`).


* `-w, --width <INT>`: Line folding/wrapping width for output sequence lines.



---

### 3. FASTQ & K-mer Analysis

#### `fq-stats`

Generate read count and quality statistics for FASTQ files.

```bash
myrs fq-stats -f reads.fastq.gz -t 4

```

* `-f, --fname <FILE>`: Input FASTQ file path (supports `.gz`).


* `-t, --threads <INT>`: Number of threads (default: `2`).



---

#### `dump-kmc`

Dump KMC database k-mers and their counts to text format (`kmer\tcount`).

```bash
myrs dump-kmc -p kmc_db_prefix -o output_kmers.txt --in-memory

```

* `-p, --kmc <PREFIX>`: Base prefix name of KMC database (e.g., `db` for `db.kmc_pre` / `db.kmc_suf`).


* `-o, --output <FILE>`: Output text file path (defaults to stdout if omitted).


* `-i, --in-memory`: Pre-load suffix buffers completely into memory instead of memory mapping (default: `false`).



---

#### `fa-kmer-cov` (not yet implemented completely)

Evaluate k-mer coverage profile across FASTA records using a KMC database.

```bash
myrs fa-kmer-cov -r ref.fasta -k kmc_db_prefix -o coverage.txt -t 4

```

* `-r, --reference <FILE>`: Input reference FASTA file path.


* `-k, --kmc <PREFIX>`: Base prefix name of KMC database.


* `-o, --output <FILE>`: Output report file path.


* `-m, --memory`: Pre-load KMC database buffers into memory (default: `false`).


* `-t, --threads <INT>`: Number of processing threads (default: `2`).



---

### 4. Chromosome Conformation (Hi-C)

#### `hic-contact-matrix`

Construct a sparse contact matrix from a Hi-C BAM file.

```bash
myrs hic-contact-matrix -i aligned_hic.bam -b 10000 -t 4 -q 30 --max-nm 5 --min-as 100

```

* `-i, --input <FILE>`: Input Hi-C BAM file path.


* `-b, --binsize <INT>`: Bin size in base pairs (e.g., `10000` for 10 kb).


* `-t, --threads <INT>`: Number of processing threads (default: `2`).


* `-q, --min-mapq <INT>`: Minimum mapping quality threshold (default: `30`).


* `--max-nm <INT>`: Maximum allowed edit distance/mismatches (`NM` tag, default: `5`).


* `--min-as <INT>`: Minimum allowed alignment score (`AS` tag, default: `100`).

---

## GenAI Involvement
This project was developed with the assistance of GenAI tools (Copilot), which helped only with drafting documentation and refactoring (along with cargo fmt). The final content has been reviewed and verified.

## Contact
For questions, issues, or feature requests, please open an issue on the [GitHub repository](https://github.com/sivasubramanics/myrs/issues).