use std::cmp::Ordering;
use std::fmt;

/// Bit-packed K-mer supporting length K up to 128 (256 bits total).
/// Uses 4 x u64 words: data[0] is MSB (5' end), data[3] is LSB (3' end).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct PackedKmer128 {
    data: [u64; 4],
    k: u8,
}

impl PackedKmer128 {
    #[inline(always)]
    pub fn new(k: usize) -> Self {
        assert!(k <= 128, "PackedKmer128 supports k <= 128");
        Self {
            data: [0, 0, 0, 0],
            k: k as u8,
        }
    }

    #[inline(always)]
    pub fn k(&self) -> usize {
        self.k as usize
    }

    #[inline(always)]
    pub fn data(&self) -> &[u64; 4] {
        &self.data
    }

    /// Appends a base byte (A, C, G, T) to the right (3' end) in O(1) time
    /// by shifting all 4 u64 words left by 2 bits.
    #[inline(always)]
    pub fn push_base(&mut self, base_byte: u8) {
        let bits = match base_byte {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00,
        };

        // Shift 256-bit word left by 2 bits
        self.data[0] = (self.data[0] << 2) | (self.data[1] >> 62);
        self.data[1] = (self.data[1] << 2) | (self.data[2] >> 62);
        self.data[2] = (self.data[2] << 2) | (self.data[3] >> 62);
        self.data[3] = (self.data[3] << 2) | (bits as u64);

        self.apply_mask();
    }

    /// Clears any bits beyond the 2*k capacity limit.
    #[inline(always)]
    fn apply_mask(&mut self) {
        let total_bits = (self.k as usize) * 2;

        // Mask top words if total bits < word boundaries
        if total_bits <= 64 {
            self.data[0] = 0;
            self.data[1] = 0;
            self.data[2] = 0;
            self.data[3] &= if total_bits == 64 { u64::MAX } else { (1u64 << total_bits) - 1 };
        } else if total_bits <= 128 {
            self.data[0] = 0;
            self.data[1] = 0;
            let b2 = total_bits - 64;
            self.data[2] &= if b2 == 64 { u64::MAX } else { (1u64 << b2) - 1 };
        } else if total_bits <= 192 {
            self.data[0] = 0;
            let b1 = total_bits - 128;
            self.data[1] &= if b1 == 64 { u64::MAX } else { (1u64 << b1) - 1 };
        } else if total_bits < 256 {
            let b0 = total_bits - 192;
            self.data[0] &= (1u64 << b0) - 1;
        }
    }

    /// Creates a K-mer from ASCII sequence bytes.
    pub fn from_ascii(sequence: &[u8]) -> Option<Self> {
        let k = sequence.len();
        if k > 128 {
            return None;
        }

        let mut kmer = Self::new(k);
        for &b in sequence {
            match b {
                b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => kmer.push_base(b),
                _ => return None,
            }
        }
        Some(kmer)
    }

    /// Reconstructs ASCII String representation.
    pub fn to_string(&self) -> String {
        const LOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let k = self.k as usize;
        let mut bytes = vec![b'A'; k];

        for i in 0..k {
            let bit_pos = (k - 1 - i) * 2;
            let word_idx = 3 - (bit_pos / 64);
            let shift = bit_pos % 64;
            let base_bits = ((self.data[word_idx] >> shift) & 0b11) as usize;
            bytes[i] = LOOKUP[base_bits];
        }

        String::from_utf8(bytes).unwrap()
    }
}

impl PartialOrd for PackedKmer128 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PackedKmer128 {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.data.cmp(&other.data)
    }
}

impl fmt::Display for PackedKmer128 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

pub struct BitPackedKmer128Iterator<'a> {
    sequence: &'a [u8],
    k: usize,
    cursor: usize,
    current_kmer: PackedKmer128,
    valid_bases: usize,
}

impl<'a> BitPackedKmer128Iterator<'a> {
    pub fn new(sequence: &'a [u8], k: usize) -> Self {
        assert!(k <= 128, "PackedKmer128Iterator supports k <= 128");
        Self {
            sequence,
            k,
            cursor: 0,
            current_kmer: PackedKmer128::new(k),
            valid_bases: 0,
        }
    }
}

impl<'a> Iterator for BitPackedKmer128Iterator<'a> {
    type Item = PackedKmer128;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.cursor < self.sequence.len() {
            let byte = self.sequence[self.cursor];
            self.cursor += 1;

            match byte {
                b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {
                    self.current_kmer.push_base(byte);
                    self.valid_bases += 1;

                    if self.valid_bases >= self.k {
                        self.valid_bases = self.k;
                        return Some(self.current_kmer);
                    }
                }
                _ => {
                    // Non-ACGT base reset
                    self.valid_bases = 0;
                    self.current_kmer = PackedKmer128::new(self.k);
                }
            }
        }
        None
    }
}