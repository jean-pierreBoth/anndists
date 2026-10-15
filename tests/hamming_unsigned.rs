use anndists::dist::{DistHamming, Distance};

macro_rules! unsigned_hamming_tests {
    ($module:ident, $unsigned:ty) => {
        mod $module {
            use super::*;

            fn check(a: &[$unsigned], b: &[$unsigned]) {
                let mismatches = a.iter().zip(b).filter(|(a, b)| a != b).count();
                let want = mismatches as f32 / a.len() as f32;
                let got = DistHamming.eval(a, b);
                if want.is_nan() {
                    assert!(got.is_nan());
                } else {
                    assert_eq!(got.to_bits(), want.to_bits(), "length {}", a.len());
                }
            }

            #[test]
            fn vector_boundaries_and_unaligned_slices() {
                for len in 0..=257 {
                    let a: Vec<$unsigned> = (0..len + 1)
                        .map(|i| {
                            (i as u64)
                                .wrapping_mul(0x9e3779b97f4a7c15)
                                .rotate_left((i % 64) as u32) as $unsigned
                        })
                        .collect();
                    let mut b = vec![0; len + 2];
                    b[2..].copy_from_slice(&a[1..]);
                    check(&a[1..], &b[2..]);
                    for (i, value) in b[2..].iter_mut().enumerate() {
                        if i % 3 != 0 {
                            *value ^= (1 as $unsigned) << (i % <$unsigned>::BITS as usize);
                        }
                    }
                    check(&a[1..], &b[2..]);
                    check(&b[2..], &a[1..]);
                }
            }

            #[test]
            fn every_hash_bit_counts_as_one_mismatch() {
                let len = 32 * <$unsigned>::BITS as usize + 3;
                let a = vec![<$unsigned>::MAX; len];
                let b: Vec<$unsigned> = (0..len)
                    .map(|i| {
                        <$unsigned>::MAX ^ ((1 as $unsigned) << (i % <$unsigned>::BITS as usize))
                    })
                    .collect();
                check(&a, &a);
                check(&a, &b);
                assert_eq!(DistHamming.eval(&a, &b), 1.0);
            }

            #[test]
            fn long_sketches_do_not_overflow_lane_counters() {
                let len = 32 * (i16::MAX as usize + 1) + 1;
                let a = vec![0; len];
                let mut b = vec![<$unsigned>::MAX; len];
                check(&a, &b);
                assert_eq!(DistHamming.eval(&a, &b), 1.0);
                for i in (0..len).step_by(3) {
                    b[i] = 0;
                }
                check(&a, &b);
            }

            #[test]
            fn empty_sketches_return_nan() {
                let empty: &[$unsigned] = &[];
                check(empty, empty);
            }

            #[test]
            #[should_panic]
            fn unequal_lengths_are_rejected() {
                DistHamming.eval(&[0 as $unsigned], &[]);
            }
        }
    };
}

unsigned_hamming_tests!(u16_sketches, u16);
unsigned_hamming_tests!(u32_sketches, u32);
unsigned_hamming_tests!(u64_sketches, u64);
