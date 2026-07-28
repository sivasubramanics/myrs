use crate::utils::defaults::COMPLEMENT;
use std::borrow::Cow;
use std::hash::Hash;

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

    /// Generates the reverse complement sequence into a new byte vector.
    pub fn revcomp(&self) -> Vec<u8> {
        let mut rev = Vec::with_capacity(self.sequence.len());
        for &byte in self.sequence.iter().rev() {
            rev.push(unsafe { *COMPLEMENT.get_unchecked(byte as usize) });
        }
        rev
    }

    /// Returns `true` if forward strand is lexicographically smaller or equal
    /// to its reverse complement.
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

    /// Zero-allocation canonical view accessor.
    /// Returns Cow::Borrowed if already canonical, or Cow::Owned if reversed.
    #[inline]
    pub fn canonical(&self) -> Cow<'a, [u8]> {
        if self.is_canonical() {
            Cow::Borrowed(self.sequence)
        } else {
            Cow::Owned(self.revcomp())
        }
    }
}

/// Zero-allocation sliding-window iterator over sequence K-mers.
/// Skips non-ACGT bases without losing continuous sliding windows.
pub struct KmerIterator<'a> {
    sequence: &'a [u8],
    k: usize,
    offset: usize,
    canonical: bool,
}

impl<'a> KmerIterator<'a> {
    #[inline]
    pub fn new(sequence: &'a [u8], k: usize, canonical: bool) -> Self {
        Self {
            sequence,
            k,
            offset: 0,
            canonical,
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = Cow<'a, [u8]>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.offset + self.k <= self.sequence.len() {
            let start = self.offset;
            let window = &self.sequence[start..start + self.k];

            // 1. Scan window for non-ACGT characters (e.g. 'N')
            let mut invalid_idx = None;
            for (idx, &byte) in window.iter().enumerate() {
                match byte {
                    b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {}
                    _ => {
                        invalid_idx = Some(idx);
                        break;
                    }
                }
            }

            // 2. If invalid character found, jump iterator past it
            if let Some(bad_pos) = invalid_idx {
                self.offset = start + bad_pos + 1;
                continue;
            }

            // Advance offset for next iteration
            self.offset += 1;

            // 3. Return canonical or raw slice based on `canonical` flag
            let kmer = Kmer::new(window);
            if self.canonical {
                return Some(kmer.canonical());
            } else {
                return Some(Cow::Borrowed(window));
            }
        }
        None
    }
}