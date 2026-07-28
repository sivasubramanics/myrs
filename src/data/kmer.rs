use crate::utils::defaults::COMPLEMENT;

#[derive(Debug, Clone, Copy)]
pub struct Kmer<'a> {
    sequence: &'a [u8],
}

impl<'a> Kmer<'a> {
    #[inline]
    pub fn new(sequence: &'a [u8]) -> Self {
        Self { sequence }
    }

    #[inline]
    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Generates the reverse complement sequence as a byte vector.
    pub fn revcomp(&self) -> Vec<u8> {
        let mut rev = Vec::with_capacity(self.sequence.len());
        for &byte in self.sequence.iter().rev() {
            rev.push(unsafe { *COMPLEMENT.get_unchecked(byte as usize) });
        }
        rev
    }

    /// Returns `true` if the forward sequence is lexicographically smaller
    /// or equal to its reverse complement.
    pub fn is_canonical(&self) -> bool {
        let len = self.sequence.len();
        for i in 0..len {
            let fwd = self.sequence[i].to_ascii_uppercase();
            let rev = COMPLEMENT[self.sequence[len - 1 - i] as usize];
            if fwd < rev {
                return true;
            } else if fwd > rev {
                return false;
            }
        }
        true
    }

    /// Returns the canonical sequence (lexicographically smallest representation
    /// between forward and reverse complement).
    pub fn canonical_bytes(&self) -> Vec<u8> {
        if self.is_canonical() {
            self.sequence.to_vec()
        } else {
            self.revcomp()
        }
    }
}

/// Zero-allocation sliding-window iterator over sequence K-mers.
pub struct KmerIterator<'a> {
    sequence: &'a [u8],
    k: usize,
    offset: usize,
    canonical_only: bool,
}

impl<'a> KmerIterator<'a> {
    #[inline]
    pub fn new(sequence: &'a [u8], k: usize, canonical_only: bool) -> Self {
        Self {
            sequence,
            k,
            offset: 0,
            canonical_only,
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = Kmer<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.offset + self.k <= self.sequence.len() {
            let window = &self.sequence[self.offset..self.offset + self.k];
            self.offset += 1;
            let kmer = Kmer::new(window);

            if self.canonical_only && !kmer.is_canonical() {
                continue; // Skip non-canonical k-mers if filtering
            }
            return Some(kmer);
        }
        None
    }
}