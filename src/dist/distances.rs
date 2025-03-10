//! Some standard distances as L1, L2, Cosine, Jaccard, Hamming
//! and a structure to enable the user to implement its own distances.
//! For the heavily used case (f32) we provide simd avx2 and std::simd implementations.

#[cfg(feature = "stdsimd")]
use super::distsimd::*;

#[cfg(feature = "simdeez_f")]
use super::disteez::*;

#[cfg(feature = "half")]
use half::f16;
use num_traits::Float;

/// The trait describing distance.
/// For example for the L1 distance
///
/// pub struct DistL1;
///
/// implement Distance<f32> for DistL1 {
/// }
///
///
/// The L1 and Cosine distance are implemented for u16, i32, i64, f32, f64
///
///
use std::os::raw::*;

#[allow(unused)]
enum DistKind {
    DistL1(String),
    DistL2(String),
    /// This is the same as Cosine dist but all data L2-normalized to 1.
    DistDot(String),
    DistCosine(String),
    DistHamming(String),
    DistJaccard(String),
    DistHellinger(String),
    DistJeffreys(String),
    DistJensenShannon(String),
    /// To store a distance defined by a C pointer function
    DistCFnPtr,
    /// Distance defined by a closure
    DistFn,
    /// Distance defined by a fn Rust pointer
    DistPtr,
    DistLevenshtein(String),
    /// used only with reloading only graph data from a previous dump
    DistNoDist(String),
}

/// This is the basic Trait describing a distance. The structure Hnsw can be instantiated by anything
/// satisfying this Trait. The crate provides implmentations for L1, L2 , Cosine, Jaccard, Hamming.
/// For other distances implement the trait possibly with the newtype pattern
pub trait Distance<T: Send + Sync> {
    fn eval(&self, va: &[T], vb: &[T]) -> f32;
}

/// Special forbidden computation distance. It is associated to a unit NoData structure
/// This is a special structure used when we want to only reload the graph from a previous computation
/// possibly from an foreign language (and we do not have access to the original type of data from the foreign language).
#[derive(Default, Copy, Clone)]
pub struct NoDist;

impl<T: Send + Sync> Distance<T> for NoDist {
    fn eval(&self, _va: &[T], _vb: &[T]) -> f32 {
        log::error!("panic error : cannot call eval on NoDist");
        panic!("cannot call distance with NoDist");
    }
} // end impl block for NoDist

/// L1 distance : implemented for i32, f64, i64, u32 , u16 , u8 and with Simd avx2 for f32
#[derive(Default, Copy, Clone)]
pub struct DistL1;

macro_rules! implementL1Distance (
    ($ty:ty) => (

    impl Distance<$ty> for DistL1  {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
            assert_eq!(va.len(), vb.len());
        // RUSTFLAGS = "-C opt-level=3 -C target-cpu=native"
            va.iter().zip(vb.iter()).map(|t| (*t.0 as f32- *t.1 as f32).abs()).sum()
        } // end of compute
    } // end of impl block
    )  // end of pattern matching
);

implementL1Distance!(i32);
implementL1Distance!(f64);
implementL1Distance!(i64);
implementL1Distance!(u32);
implementL1Distance!(u16);
implementL1Distance!(u8);

impl Distance<f32> for DistL1 {
    fn eval(&self, va: &[f32], vb: &[f32]) -> f32 {
        //
        cfg_if::cfg_if! {
        if #[cfg(feature = "simdeez_f")] {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))] {
                if is_x86_feature_detected!("avx2") {
                    return unsafe {distance_l1_f32_avx2(va,vb)};
                }
                else {
                    assert_eq!(va.len(), vb.len());
                    va.iter().zip(vb.iter()).map(|t| (*t.0 as f32- *t.1 as f32).abs()).sum()
                }
            }
        }
        else if #[cfg(feature = "stdsimd")] {
            distance_l1_f32_simd(va,vb)
        }
        else {
            va.iter().zip(vb.iter()).map(|t| (*t.0 - *t.1 ).abs()).sum()
        }
        } // end cfg_if
    } // end of eval
} // end impl Distance<f32> for DistL1

// L1 for f16
#[cfg(feature = "half")]
impl Distance<f16> for DistL1 {
    fn eval(&self, va: &[f16], vb: &[f16]) -> f32 {
        assert_eq!(va.len(), vb.len());
        va.iter()
            .zip(vb.iter())
            .map(|(a, b)| Float::abs(a - b))
            .sum::<f16>()
            .to_f32()
    }
}

//========================================================================

/// L2 distance : implemented for i32, f64, i64, u32 , u16 , u8 and with Simd avx2 for f32
#[derive(Default, Copy, Clone)]
pub struct DistL2;

macro_rules! implementL2Distance (
    ($ty:ty) => (

    impl Distance<$ty> for DistL2  {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
            assert_eq!(va.len(), vb.len());
            let norm : f32 = va.iter().zip(vb.iter()).map(|t| (*t.0 as f32- *t.1 as f32) * (*t.0 as f32- *t.1 as f32)).sum();
            norm.sqrt()
        } // end of compute
    } // end of impl block
    )  // end of pattern matching
);

//implementL2Distance!(f32);
implementL2Distance!(i32);
implementL2Distance!(f64);
implementL2Distance!(i64);
implementL2Distance!(u32);
implementL2Distance!(u16);
implementL2Distance!(u8);

#[allow(unused)]
// base scalar l2 for f32
fn scalar_l2_f32(va: &[f32], vb: &[f32]) -> f32 {
    let norm: f32 = va
        .iter()
        .zip(vb.iter())
        .map(|t| (*t.0 - *t.1) * (*t.0 - *t.1))
        .sum();
    assert!(norm >= 0.);
    norm.sqrt()
}

impl Distance<f32> for DistL2 {
    fn eval(&self, va: &[f32], vb: &[f32]) -> f32 {
        //
        cfg_if::cfg_if! {
            if #[cfg(feature = "simdeez_f")] {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                if is_x86_feature_detected!("avx2") {
                    return unsafe { distance_l2_f32_avx2(va, vb) };
                }
                else {
                    return scalar_l2_f32(&va, &vb);
                }
            }
            } else if #[cfg(feature = "stdsimd")] {
                return distance_l2_f32_simd(va, vb);
            }
            else {
                scalar_l2_f32(va, vb)
            }
        }
    } // end of eval
} // end impl Distance<f32> for DistL2

//=========================================================================

/// Cosine distance : implemented for f16, f32, f64, i64, i32 , u16
#[derive(Default, Copy, Clone)]
pub struct DistCosine;

macro_rules! implementCosDistance(
    ($ty:ty) => (
     impl Distance<$ty> for DistCosine  {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
            assert_eq!(va.len(), vb.len());
            //
            let dist:f32;
            let zero:f64 = 0.;
            // to // by rayon
            let res = va.iter().zip(vb.iter()).map(|t| ((*t.0 * *t.1) as f64, (*t.0 * *t.0) as f64, (*t.1 * *t.1) as f64)).
                fold((0., 0., 0.), |acc , t| (acc.0 + t.0, acc.1 + t.1, acc.2 + t.2));
            //
            if res.1 > zero && res.2 > zero {
                let dist_unchecked = 1. - res.0 / (res.1 * res.2).sqrt();
                assert!(dist_unchecked >= - 0.00002);
                dist = dist_unchecked.max(0.) as f32;
            }
            else {
                dist = 0.;
            }
            //
            return dist;
        } // end of function
     } // end of impl block
    ) // end of matching
);

implementCosDistance!(f32);
implementCosDistance!(f64);
implementCosDistance!(i64);
implementCosDistance!(i32);
implementCosDistance!(u16);

// Cosine for f16
#[cfg(feature = "half")]
impl Distance<f16> for DistCosine {
    fn eval(&self, va: &[f16], vb: &[f16]) -> f32 {
        assert_eq!(va.len(), vb.len());
        let res = va
            .iter()
            .zip(vb.iter())
            .map(|(a, b)| ((a * b), (a * a), (b * b)))
            .fold((0., 0., 0.), |acc, (ab, aa, bb)| {
                (
                    acc.0 + ab.to_f32(),
                    acc.1 + aa.to_f32(),
                    acc.2 + bb.to_f32(),
                )
            });
        if res.1 > 0. && res.2 > 0. {
            let dist_unchecked = 1. - res.0 / (res.1 * res.2).sqrt();
            assert!(dist_unchecked >= -0.00002);
            dist_unchecked.max(0.)
        } else {
            0.
        }
    }
}

//=========================================================================

/// This is essentially the Cosine distance but we suppose
/// all vectors (graph construction and request vectors have been l2 normalized to unity
/// BEFORE INSERTING in  HNSW!.   
/// No control is made, so it is the user responsability to send normalized vectors
/// everywhere in inserting and searching.
///
/// In large dimensions (hundreds) this pre-normalization spare cpu time.  
/// At low dimensions (a few ten's there is not a significant gain).  
/// This distance makes sense only for f16, f32 or f64
/// We provide for avx2 implementations for f32 that provides consequent gains
/// in large dimensions

#[derive(Default, Copy, Clone)]
pub struct DistDot;

#[allow(unused)]
macro_rules! implementDotDistance(
    ($ty:ty) => (
     impl Distance<$ty> for DistDot  {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
            assert_eq!(va.len(), vb.len());
            //
            let zero:f32 = 0f32;
            // to // by rayon
            let dot = va.iter().zip(vb.iter()).map(|t| (*t.0 * *t.1) as f32).fold(0., |acc , t| (acc + t));
            //
            assert(dot <= 1.);
            return  1. - dot;
        } // end of function
      } // end of impl block
    ) // end of matching
);

#[allow(unused)]
fn scalar_dot_f32(va: &[f32], vb: &[f32]) -> f32 {
    let dot = 1.
        - va.iter()
            .zip(vb.iter())
            .map(|t| (*t.0 * *t.1))
            .fold(0., |acc, t| (acc + t));
    assert!(dot >= 0.);
    dot
}

impl Distance<f32> for DistDot {
    fn eval(&self, va: &[f32], vb: &[f32]) -> f32 {
        //
        cfg_if::cfg_if! {
            if #[cfg(feature = "simdeez_f")] {
                #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                {
                    if is_x86_feature_detected!("avx2") {
                        return unsafe { distance_dot_f32_avx2(va, vb) };
                    } else if is_x86_feature_detected!("sse2") {
                        return unsafe { distance_dot_f32_sse2(va, vb) };
                    }
                    else {
                        return scalar_dot_f32(va, vb);
                    }
                } // end x86
            } else if #[cfg(feature = "stdsimd")] {
                return distance_dot_f32_simd_iter(va,vb);
            }
            else {
                scalar_dot_f32(va, vb)
            }
        }
    } // end of eval
}

// Dot distance for f16
#[cfg(feature = "half")]
impl Distance<f16> for DistDot {
    fn eval(&self, va: &[f16], vb: &[f16]) -> f32 {
        assert_eq!(va.len(), vb.len());
        let dot: f32 = va
            .iter()
            .zip(vb.iter())
            .map(|(a, b)| a * b)
            .sum::<f16>()
            .to_f32();
        1. - dot
    }
}

pub fn l2_normalize(va: &mut [f32]) {
    let l2norm = va.iter().map(|t| (*t * *t)).sum::<f32>().sqrt();
    if l2norm > 0. {
        for v in va {
            *v /= l2norm;
        }
    }
}

//=======================================================================================

///
/// A structure to compute Hellinger distance between probalilities.
/// Vector must be >= 0 and normalized to 1.
///   
/// The distance computation does not check that
/// and in fact simplifies the expression of distance assuming vectors are positive and L1 normalised to 1.
/// The user must enforce these conditions before  inserting otherwise results will be meaningless
/// at best or code will panic!
///
/// For f32 a simd implementation is provided if avx2 is detected.
#[derive(Default, Copy, Clone)]
pub struct DistHellinger;

// default implementation
macro_rules! implementHellingerDistance (
    ($ty:ty) => (

    impl Distance<$ty> for DistHellinger {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
            assert_eq!(va.len(), vb.len());
        // RUSTFLAGS = "-C opt-level=3 -C target-cpu=native"
        // to // by rayon
            let mut dist = va.iter().zip(vb.iter()).map(|t| ((*t.0).sqrt() * (*t.1).sqrt()) as f32).fold(0., |acc , t| (acc + t*t));
            dist = (1. - dist).sqrt();
            dist
        } // end of compute
    } // end of impl block
    )  // end of pattern matching
);

implementHellingerDistance!(f64);

impl Distance<f32> for DistHellinger {
    fn eval(&self, va: &[f32], vb: &[f32]) -> f32 {
        //
        #[cfg(feature = "simdeez_f")]
        {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                if is_x86_feature_detected!("avx2") {
                    return unsafe { distance_hellinger_f32_avx2(va, vb) };
                }
            }
        }
        let mut dist = va
            .iter()
            .zip(vb.iter())
            .map(|t| ((*t.0) * (*t.1)).sqrt())
            .fold(0., |acc, t| (acc + t));
        // if too far away from >= panic else reset!
        assert!(1. - dist >= -0.000001);
        dist = (1. - dist).max(0.).sqrt();
        dist
    } // end of eval
}

// Hellinger for f16
#[cfg(feature = "half")]
impl Distance<f16> for DistHellinger {
    fn eval(&self, va: &[f16], vb: &[f16]) -> f32 {
        assert_eq!(va.len(), vb.len());
        let mut dist = va
            .iter()
            .zip(vb.iter())
            .map(|(a, b)| a.sqrt() * b.sqrt())
            .sum::<f16>()
            .to_f32();
        assert!(1. - dist >= -0.000001);
        dist = (1. - dist).max(0.).sqrt();
        dist
    }
}

//=======================================================================================

///
/// A structure to compute Jeffreys divergence between probalilities.
/// If p and q are 2 probability distributions
/// the "distance" is computed as:
///   sum (p\[i\] - q\[i\]) * ln(p\[i\]/q\[i\])
///
/// To take care of null probabilities in the formula we use  max(x\[i\],1.E-30)
/// for x = p and q in the log compuations
///   
/// Vector must be >= 0 and normalized to 1!  
/// The distance computation does not check that.
/// The user must enforce these conditions before inserting in the hnws structure,
/// otherwise results will be meaningless at best or code will panic!
///
/// For f32 a simd implementation is provided if avx2 is detected.
#[derive(Default, Copy, Clone)]
pub struct DistJeffreys;

pub const M_MIN: f32 = 1.0e-30;

// default implementation
macro_rules! implementJeffreysDistance (
    ($ty:ty) => (

    impl Distance<$ty> for DistJeffreys {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
        // RUSTFLAGS = "-C opt-level=3 -C target-cpu=native"
        let dist = va.iter().zip(vb.iter()).map(|t| (*t.0 - *t.1) * ((*t.0).max(M_MIN as f64)/ (*t.1).max(M_MIN as f64)).ln() as f64).fold(0., |acc , t| (acc + t*t));
        dist as f32
        } // end of compute
    } // end of impl block
    )  // end of pattern matching
);

implementJeffreysDistance!(f64);

impl Distance<f32> for DistJeffreys {
    fn eval(&self, va: &[f32], vb: &[f32]) -> f32 {
        //
        #[cfg(feature = "simdeez_f")]
        {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                if is_x86_feature_detected!("avx2") {
                    return unsafe { distance_jeffreys_f32_avx2(va, vb) };
                }
            }
        }
        let dist = va
            .iter()
            .zip(vb.iter())
            .map(|t| (*t.0 - *t.1) * ((*t.0).max(M_MIN) / (*t.1).max(M_MIN)).ln())
            .fold(0., |acc, t| (acc + t));
        dist
    } // end of eval
}

// Jeffreys for f16
#[cfg(feature = "half")]
impl Distance<f16> for DistJeffreys {
    fn eval(&self, va: &[f16], vb: &[f16]) -> f32 {
        assert_eq!(va.len(), vb.len());
        va.iter()
            .zip(vb.iter())
            .map(|(a, b)| {
                let (a_clamp, b_clamp) = (a.max(f16::MIN), b.max(f16::MIN));
                (a - b) * (a_clamp / b_clamp).ln()
            })
            .sum::<f16>()
            .to_f32()
    }
}

//=======================================================================================

/// Jensen-Shannon distance.  
/// It is defined as the **square root** of the  Jensen–Shannon divergence and is a metric.
/// Vector must be >= 0 and normalized to 1!
/// **The distance computation does not check that**.
#[derive(Default, Copy, Clone)]
pub struct DistJensenShannon;

macro_rules! implementDistJensenShannon (

    ($ty:ty) => (
        impl Distance<$ty> for DistJensenShannon {
            fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
                let mut dist = 0.;
                //
                assert_eq!(va.len(), vb.len());
                //
                for i in 0..va.len() {
                    let mean_ab = 0.5 * (va[i] + vb[i]);
                    if va[i] > 0. {
                        dist += va[i] * (va[i]/mean_ab).ln();
                    }
                    if vb[i] > 0. {
                        dist += vb[i] * (vb[i]/mean_ab).ln();
                    }
                }
                (0.5 * dist).sqrt() as f32
            } // end eval
        }  // end impl Distance<$ty>
    )  // end of pattern matching on ty
);

implementDistJensenShannon!(f64);
implementDistJensenShannon!(f32);

/// Jensen-Shannon for f16
#[cfg(feature = "half")]
impl Distance<f16> for DistJensenShannon {
    fn eval(&self, va: &[f16], vb: &[f16]) -> f32 {
        assert_eq!(va.len(), vb.len());
        let mut dist = 0_f32;
        for i in 0..va.len() {
            let (af, bf) = (va[i].to_f32(), vb[i].to_f32());
            let mean_ab = 0.5 * (af + bf);
            if af > 0. {
                dist += af * (af / mean_ab).ln();
            }
            if bf > 0. {
                dist += bf * (bf / mean_ab).ln();
            }
        }
        (0.5 * dist).sqrt()
    }
}

//=======================================================================================

/// Hamming distance. Implemented for u8, u16, u32, i32 and i16
/// The distance returned is normalized by length of slices, so it is between 0. and 1.  
///
/// A special implementation for f64 is made but exclusively dedicated to SuperMinHash usage in crate [probminhash](https://crates.io/crates/probminhash).  
/// It could be made generic with the PartialEq implementation for f64 and f32 in unsable source of Rust
#[derive(Default, Copy, Clone)]
pub struct DistHamming;

macro_rules! implementHammingDistance (
    ($ty:ty) => (

    impl Distance<$ty> for DistHamming  {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
        // RUSTFLAGS = "-C opt-level=3 -C target-cpu=native"
            assert_eq!(va.len(), vb.len());
            let norm : f32 = va.iter().zip(vb.iter()).filter(|t| t.0 != t.1).count() as f32;
            norm / va.len() as f32
        } // end of compute
    } // end of impl block
    )  // end of pattern matching
);

impl Distance<i32> for DistHamming {
    fn eval(&self, va: &[i32], vb: &[i32]) -> f32 {
        //
        #[cfg(feature = "simdeez_f")]
        {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                if is_x86_feature_detected!("avx2") {
                    return unsafe { distance_hamming_i32_avx2(va, vb) };
                }
            }
        }
        assert_eq!(va.len(), vb.len());
        let dist: f32 = va.iter().zip(vb.iter()).filter(|t| t.0 != t.1).count() as f32;
        dist / va.len() as f32
    } // end of eval
} // end implementation Distance<i32>

/// This implementation is dedicated to SuperMinHash algorithm in crate [probminhash](https://crates.io/crates/probminhash).  
/// Could be made generic with unstable source as there is implementation of PartialEq for f64
impl Distance<f64> for DistHamming {
    fn eval(&self, va: &[f64], vb: &[f64]) -> f32 {
        /*   Tests show that it is slower than basic method!!!
        #[cfg(feature = "simdeez_f")] {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))] {
                if is_x86_feature_detected!("avx2") {
                    log::trace!("calling distance_hamming_f64_avx2");
                    return unsafe { distance_hamming_f64_avx2(va,vb) };
                }
            }
        }
        */
        //
        assert_eq!(va.len(), vb.len());
        let dist: usize = va.iter().zip(vb.iter()).filter(|t| t.0 != t.1).count();
        (dist as f64 / va.len() as f64) as f32
    } // end of eval
} // end implementation Distance<f64>

//

/// This implementation is dedicated to SuperMinHash algorithm in crate [probminhash](https://crates.io/crates/probminhash).  
/// Could be made generic with unstable source as there is implementation of PartialEq for f32
impl Distance<f32> for DistHamming {
    fn eval(&self, va: &[f32], vb: &[f32]) -> f32 {
        cfg_if::cfg_if! {
            if #[cfg(feature = "stdsimd")] {
                return distance_jaccard_f32_16_simd(va,vb);
            }
            else {
                assert_eq!(va.len(), vb.len());
                let dist : usize = va.iter().zip(vb.iter()).filter(|t| t.0 != t.1).count();
                (dist as f64 / va.len() as f64) as f32
            }
        }
    } // end of eval
} // end implementation Distance<f32>

//

#[cfg(feature = "stdsimd")]
impl Distance<u32> for DistHamming {
    fn eval(&self, va: &[u32], vb: &[u32]) -> f32 {
        //
        return distance_jaccard_u32_16_simd(va, vb);
    } // end of eval
} // end implementation Distance<u32>

//

#[cfg(feature = "stdsimd")]
impl Distance<u64> for DistHamming {
    fn eval(&self, va: &[u64], vb: &[u64]) -> f32 {
        return distance_jaccard_u64_8_simd(va, vb);
    } // end of eval
} // end implementation Distance<u64>

// i32 is implmeented by simd
implementHammingDistance!(u8);
implementHammingDistance!(u16);

#[cfg(not(feature = "stdsimd"))]
implementHammingDistance!(u32);

#[cfg(not(feature = "stdsimd"))]
implementHammingDistance!(u64);

implementHammingDistance!(i16);

//====================================================================================
//   Jaccard Distance

/// Jaccard distance. Implemented for u8, u16 , u32.
#[derive(Default, Copy, Clone)]
pub struct DistJaccard;

// contruct a 2-uple accumulator that has sum of max in first component , and sum of min in 2 component
// stay in integer as long as possible
// Note : summing u32 coming from hash values can overflow! We must go up to u64 for additions!
macro_rules! implementJaccardDistance (
    ($ty:ty) => (

    impl Distance<$ty> for DistJaccard  {
        fn eval(&self, va:&[$ty], vb: &[$ty]) -> f32 {
            let (max,min) : (u64, u64) = va.iter().zip(vb.iter()).fold((0u64,0u64), |acc, t| if t.0 > t.1 {
                                (acc.0 + *t.0 as u64, acc.1 + *t.1 as u64) }
                        else {
                                (acc.0 + *t.1 as u64 , acc.1 + *t.0 as u64)
                             }
            );
            if max > 0 {
                let dist = 1. - (min  as f64)/ (max as f64);
                assert!(dist >= 0.);
                dist as f32
            }
            else {
                0.
            }
        } // end of compute
    } // end of impl block
    )  // end of pattern matching
);

implementJaccardDistance!(u8);
implementJaccardDistance!(u16);
implementJaccardDistance!(u32);

// ==========================================================================================

/// Levenshtein distance. Implemented for u16
#[derive(Default, Copy, Clone)]
pub struct DistLevenshtein;

impl Distance<u16> for DistLevenshtein {
    fn eval(&self, a: &[u16], b: &[u16]) -> f32 {
        let len_a = a.len();
        let len_b = b.len();
        if len_a < len_b {
            return self.eval(b, a);
        }
        // handle special case of 0 length
        if len_a == 0 {
            return len_b as f32;
        } else if len_b == 0 {
            return len_a as f32;
        }

        let len_b = len_b + 1;

        let mut pre;
        let mut tmp;
        let mut cur: Vec<usize> = (0..len_b).collect();

        // calculate edit distance
        for (i, ca) in a.iter().enumerate() {
            // get first column for this row
            pre = cur[0];
            cur[0] = i + 1;
            for (j, cb) in b.iter().enumerate() {
                tmp = cur[j + 1];
                cur[j + 1] = std::cmp::min(
                    // deletion
                    tmp + 1,
                    std::cmp::min(
                        // insertion
                        cur[j] + 1,
                        // match or substitution
                        pre + if ca == cb { 0 } else { 1 },
                    ),
                );
                pre = tmp;
            }
        }
        cur[len_b - 1] as f32
    }
}

//=======================================================================================
//   Case of function pointers (cover Trait Fn , FnOnce ...)
// The book (Function item types):  " There is a coercion from function items to function pointers with the same signature  "
// The book (Call trait and coercions): "Non capturing closures can be coerced to function pointers with the same signature"

/// This type is for function with a C-API
/// Distances can be computed by such a function. It
/// takes as arguments the two (C, rust, julia) pointers to primitive type vectos and length
/// passed as a unsignedlonlong (64 bits) which is called c_ulonglong in Rust and Culonglong in Julia
///
type DistCFnPtr<T> = extern "C" fn(*const T, *const T, len: c_ulonglong) -> f32;

/// A structure to implement Distance Api for type DistCFnPtr\<T\>,
/// i.e distance provided by a C function pointer.  
/// It must be noted that this can be used in Julia via the macro @cfunction
/// to define interactiveley a distance function , compile it on the fly and sent it
/// to Rust via the init_hnsw_{f32, i32, u16, u32, u8} function
/// defined in libext
///
pub struct DistCFFI<T: Copy + Clone + Sized + Send + Sync> {
    dist_function: DistCFnPtr<T>,
}

impl<T: Copy + Clone + Sized + Send + Sync> DistCFFI<T> {
    pub fn new(f: DistCFnPtr<T>) -> Self {
        DistCFFI { dist_function: f }
    }
}

impl<T: Copy + Clone + Sized + Send + Sync> Distance<T> for DistCFFI<T> {
    fn eval(&self, va: &[T], vb: &[T]) -> f32 {
        // get pointers
        let len = va.len();
        let ptr_a = va.as_ptr();
        let ptr_b = vb.as_ptr();
        let dist = (self.dist_function)(ptr_a, ptr_b, len as c_ulonglong);
        log::trace!(
            "DistCFFI dist_function_ptr {:?} returning {:?} ",
            self.dist_function,
            dist
        );
        dist
    } // end of compute
} // end of impl block

//========================================================================================================

/// This structure is to let user define their own distance with closures.
#[allow(clippy::type_complexity)]
pub struct DistFn<T: Copy + Clone + Sized + Send + Sync> {
    dist_function: Box<dyn Fn(&[T], &[T]) -> f32 + Send + Sync>,
}

#[allow(clippy::type_complexity)]
impl<T: Copy + Clone + Sized + Send + Sync> DistFn<T> {
    /// construction of a DistFn
    pub fn new(f: Box<dyn Fn(&[T], &[T]) -> f32 + Send + Sync>) -> Self {
        DistFn { dist_function: f }
    }
}

impl<T: Copy + Clone + Sized + Send + Sync> Distance<T> for DistFn<T> {
    fn eval(&self, va: &[T], vb: &[T]) -> f32 {
        (self.dist_function)(va, vb)
    }
}

//=======================================================================================

/// This structure uses a Rust function pointer to define the distance.
/// For commodity it can build upon a fonction returning a f64.
/// Beware that if F is f64, the distance converted to f32 can overflow!

#[derive(Copy, Clone)]
pub struct DistPtr<T: Copy + Clone + Sized + Send + Sync, F: Float> {
    dist_function: fn(&[T], &[T]) -> F,
}

impl<T: Copy + Clone + Sized + Send + Sync, F: Float> DistPtr<T, F> {
    /// construction of a DistPtr
    pub fn new(f: fn(&[T], &[T]) -> F) -> Self {
        DistPtr { dist_function: f }
    }
}

/// beware that if F is f64, the distance converted to f32 can overflow!
impl<T: Copy + Clone + Sized + Send + Sync, F: Float> Distance<T> for DistPtr<T, F> {
    fn eval(&self, va: &[T], vb: &[T]) -> f32 {
        (self.dist_function)(va, vb).to_f32().unwrap()
    }
}

//=======================================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn init_log() -> u64 {
        let mut builder = env_logger::Builder::from_default_env();
        let _ = builder.is_test(true).try_init();
        println!("\n ************** initializing logger *****************\n");
        1
    }

    #[test]
    fn test_access_to_dist_l1() {
        let distl1 = DistL1;
        //
        let v1: Vec<i32> = vec![1, 2, 3];
        let v2: Vec<i32> = vec![2, 2, 3];

        let d1 = Distance::eval(&distl1, &v1, &v2);
        assert_eq!(d1, 1_f32);

        let v3: Vec<f32> = vec![1., 2., 3.];
        let v4: Vec<f32> = vec![2., 2., 3.];
        let d2 = distl1.eval(&v3, &v4);
        assert_eq!(d2, 1_f32);
    }

    #[test]
    fn have_avx2() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                println!("I have avx2");
            } else {
                println!(" ************ I DO NOT  have avx2  ***************");
            }
        }
    } // end if

    #[test]
    fn have_avx512f() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx512f") {
                println!("I have avx512f");
            } else {
                println!(" ************ I DO NOT  have avx512f  ***************");
            }
        } // end of have_avx512f
    }

    #[test]
    fn have_sse2() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("sse2") {
                println!("I have sse2");
            } else {
                println!(" ************ I DO NOT  have SSE2  ***************");
            }
        }
    } // end of have_sse2

    #[test]
    fn test_access_to_dist_cos() {
        let distcos = DistCosine;
        //
        let v1: Vec<i32> = vec![1, -1, 1];
        let v2: Vec<i32> = vec![2, 1, -1];

        let d1 = Distance::eval(&distcos, &v1, &v2);
        assert_eq!(d1, 1_f32);
        //
        let v1: Vec<f32> = vec![1.234, -1.678, 1.367];
        let v2: Vec<f32> = vec![4.234, -6.678, 10.367];
        let d1 = Distance::eval(&distcos, &v1, &v2);

        let mut normv1 = 0.;
        let mut normv2 = 0.;
        let mut prod = 0.;
        for i in 0..v1.len() {
            prod += v1[i] * v2[i];
            normv1 += v1[i] * v1[i];
            normv2 += v2[i] * v2[i];
        }
        let dcos = 1. - prod / (normv1 * normv2).sqrt();
        println!("dist cos avec macro = {:?} ,  avec for {:?}", d1, dcos);
    }

    #[test]
    fn test_dot_distances() {
        let mut v1: Vec<f32> = vec![1.234, -1.678, 1.367];
        let mut v2: Vec<f32> = vec![4.234, -6.678, 10.367];

        let mut normv1 = 0.;
        let mut normv2 = 0.;
        let mut prod = 0.;
        for i in 0..v1.len() {
            prod += v1[i] * v2[i];
            normv1 += v1[i] * v1[i];
            normv2 += v2[i] * v2[i];
        }
        let dcos = 1. - prod / (normv1 * normv2).sqrt();
        //
        l2_normalize(&mut v1);
        l2_normalize(&mut v2);

        println!(" after normalisation v1 = {:?}", v1);
        let dot = DistDot.eval(&v1, &v2);
        println!(
            "dot  cos avec prenormalisation  = {:?} ,  avec for {:?}",
            dot, dcos
        );
    }

    #[test]
    fn test_l1() {
        init_log();
        //
        let va: Vec<f32> = vec![1.234, -1.678, 1.367, 1.234, -1.678, 1.367];
        let vb: Vec<f32> = vec![4.234, -6.678, 10.367, 1.234, -1.678, 1.367];
        //
        let dist = DistL1.eval(&va, &vb);
        let dist_check = va
            .iter()
            .zip(vb.iter())
            .map(|t| (*t.0 - *t.1).abs())
            .sum::<f32>();
        //
        log::info!(" dist : {:.5e} dist_check : {:.5e}", dist, dist_check);
        assert!((dist - dist_check).abs() / dist_check < 1.0e-5);
    } // end of test_l1

    #[test]
    fn test_jaccard_u16() {
        let v1: Vec<u16> = vec![1, 2, 1, 4, 3];
        let v2: Vec<u16> = vec![2, 2, 1, 5, 6];

        let dist = DistJaccard.eval(&v1, &v2);
        println!("dist jaccard = {:?}", dist);
        assert_eq!(dist, 1. - 11. / 16.);
    } // end of test_jaccard

    #[test]
    fn test_levenshtein() {
        let mut v1: Vec<u16> = vec![1, 2, 3, 4];
        let mut v2: Vec<u16> = vec![1, 2, 3, 3];
        let mut dist = DistLevenshtein.eval(&v1, &v2);
        println!("dist levenshtein = {:?}", dist);
        assert_eq!(dist, 1.0);
        v1 = vec![1, 2, 3, 4];
        v2 = vec![1, 2, 3, 4];
        dist = DistLevenshtein.eval(&v1, &v2);
        println!("dist levenshtein = {:?}", dist);
        assert_eq!(dist, 0.0);
        v1 = vec![1, 1, 1, 4];
        v2 = vec![1, 2, 3, 4];
        dist = DistLevenshtein.eval(&v1, &v2);
        println!("dist levenshtein = {:?}", dist);
        assert_eq!(dist, 2.0);
        v2 = vec![1, 1, 1, 4];
        v1 = vec![1, 2, 3, 4];
        dist = DistLevenshtein.eval(&v1, &v2);
        println!("dist levenshtein = {:?}", dist);
        assert_eq!(dist, 2.0);
    } // end of test_levenshtein

    extern "C" fn dist_func_float(va: *const f32, vb: *const f32, len: c_ulonglong) -> f32 {
        let mut dist: f32 = 0.;
        let sa = unsafe { std::slice::from_raw_parts(va, len as usize) };
        let sb = unsafe { std::slice::from_raw_parts(vb, len as usize) };

        for i in 0..len {
            dist += (sa[i as usize] - sb[i as usize]).abs().sqrt();
        }
        dist
    }

    #[test]
    fn test_dist_ext_float() {
        let va: Vec<f32> = vec![1., 2., 3.];
        let vb: Vec<f32> = vec![1., 2., 3.];
        println!("in test_dist_ext_float");
        let dist1 = dist_func_float(va.as_ptr(), vb.as_ptr(), va.len() as c_ulonglong);
        println!("test_dist_ext_float computed : {:?}", dist1);

        let mydist = DistCFFI::<f32>::new(dist_func_float);

        let dist2 = mydist.eval(&va, &vb);
        assert_eq!(dist1, dist2);
    } // end test_dist_ext_float

    #[test]

    fn test_my_closure() {
        //        use hnsw_rs::dist::Distance;
        let weight = [0.1, 0.8, 0.1];
        let my_fn = move |va: &[f32], vb: &[f32]| -> f32 {
            // should check that we work with same size for va, vb, and weight...
            let mut dist: f32 = 0.;
            for i in 0..va.len() {
                dist += weight[i] * (va[i] - vb[i]).abs();
            }
            dist
        };
        let my_boxed_f = Box::new(my_fn);
        let my_boxed_dist = DistFn::<f32>::new(my_boxed_f);
        let va: Vec<f32> = vec![1., 2., 3.];
        let vb: Vec<f32> = vec![2., 2., 4.];
        let dist = my_boxed_dist.eval(&va, &vb);
        println!("test_my_closure computed : {:?}", dist);
        // try allocation Hnsw
        //        let _hnsw = Hnsw::<f32, hnsw_rs::dist::DistFn<f32>>::new(10, 3, 100, 16, my_boxed_dist);
        //
        assert_eq!(dist, 0.2);
    } // end of test_my_closure

    #[test]
    fn test_hellinger() {
        let length = 9;
        let mut p_data = Vec::with_capacity(length);
        let mut q_data = Vec::with_capacity(length);
        for _ in 0..length {
            p_data.push(1. / length as f32);
            q_data.push(1. / length as f32);
        }
        p_data[0] -= 1. / (2 * length) as f32;
        p_data[1] += 1. / (2 * length) as f32;
        //
        let dist = DistHellinger.eval(&p_data, &q_data);

        let dist_exact_fn = |n: usize| -> f32 {
            let d1 = (4. - (6_f32).sqrt() - (2_f32).sqrt()) / n as f32;
            d1.sqrt() / (2_f32).sqrt()
        };
        let dist_exact = dist_exact_fn(length);
        //
        log::info!("dist computed {:?} dist exact{:?} ", dist, dist_exact);
        println!("dist computed  {:?} , dist exact {:?} ", dist, dist_exact);
        //
        assert!((dist - dist_exact).abs() < 1.0e-5);
    }

    #[test]
    fn test_jeffreys() {
        // this essentially test av2 implementation for f32
        let length = 19;
        let mut p_data: Vec<f32> = Vec::with_capacity(length);
        let mut q_data: Vec<f32> = Vec::with_capacity(length);
        for _ in 0..length {
            p_data.push(1. / length as f32);
            q_data.push(1. / length as f32);
        }
        p_data[0] -= 1. / (2 * length) as f32;
        p_data[1] += 1. / (2 * length) as f32;
        q_data[10] += 1. / (2 * length) as f32;
        //
        let dist_eval = DistJeffreys.eval(&p_data, &q_data);
        let mut dist_test = 0.;
        for i in 0..length {
            dist_test +=
                (p_data[i] - q_data[i]) * (p_data[i].max(M_MIN) / q_data[i].max(M_MIN)).ln();
        }
        //
        log::info!("dist eval {:?} dist test{:?} ", dist_eval, dist_test);
        println!("dist eval  {:?} , dist test {:?} ", dist_eval, dist_test);
        assert!(dist_test >= 0.);
        assert!((dist_eval - dist_test).abs() < 1.0e-5);
    }

    #[test]
    fn test_jensenshannon() {
        init_log();
        //
        let length = 19;
        let mut p_data: Vec<f32> = Vec::with_capacity(length);
        let mut q_data: Vec<f32> = Vec::with_capacity(length);
        for _ in 0..length {
            p_data.push(1. / length as f32);
            q_data.push(1. / length as f32);
        }
        p_data[0] -= 1. / (2 * length) as f32;
        p_data[1] += 1. / (2 * length) as f32;
        q_data[10] += 1. / (2 * length) as f32;
        p_data[12] = 0.;
        q_data[12] = 0.;
        //
        let dist_eval = DistJensenShannon.eval(&p_data, &q_data);
        //
        log::info!("dist eval {:?} ", dist_eval);
        println!("dist eval  {:?} ", dist_eval);
    }

    #[allow(unused)]
    use rand::distributions::{Distribution, Uniform};

    // to be run with and without simdeez_f
    #[test]
    fn test_hamming_f64() {
        init_log();

        let size_test = 500;
        let fmax: f64 = 3.;
        let mut rng = rand::thread_rng();
        for i in 300..size_test {
            // generer 2 va et vb s des vecteurs<i32> de taille i  avec des valeurs entre -imax et + imax et controler les resultat
            let between = Uniform::<f64>::from(-fmax..fmax);
            let va: Vec<f64> = (0..i).map(|_| between.sample(&mut rng)).collect();
            let mut vb: Vec<f64> = (0..i).map(|_| between.sample(&mut rng)).collect();
            // reset half of vb to va
            vb[..(i / 2)].copy_from_slice(&va[..(i / 2)]);

            let easy_dist: u32 = va
                .iter()
                .zip(vb.iter())
                .map(|(a, b)| if a != b { 1 } else { 0 })
                .sum();
            let h_dist = DistHamming.eval(&va, &vb);
            let easy_dist = easy_dist as f32 / va.len() as f32;
            let j_exact = ((i / 2) as f32) / (i as f32);
            log::debug!(
                "test size {:?}  HammingDist {:.3e} easy : {:.3e} exact : {:.3e} ",
                i,
                h_dist,
                easy_dist,
                j_exact
            );
            if (easy_dist - h_dist).abs() > 1.0e-5 {
                println!(" jhamming = {:?} , jexact = {:?}", h_dist, easy_dist);
                log::debug!("va = {:?}", va);
                log::debug!("vb = {:?}", vb);
                std::process::exit(1);
            }
            if (j_exact - h_dist).abs() > 1. / i as f32 + 1.0E-5 {
                println!(
                    " jhamming = {:?} , jexact = {:?}, j_easy : {:?}",
                    h_dist, j_exact, easy_dist
                );
                log::debug!("va = {:?}", va);
                log::debug!("vb = {:?}", vb);
                std::process::exit(1);
            }
        }
    } // end of test_hamming_f64

    #[test]
    fn test_hamming_f32() {
        init_log();

        let size_test = 500;
        let fmax: f32 = 3.;
        let mut rng = rand::thread_rng();
        for i in 300..size_test {
            // generer 2 va et vb s des vecteurs<i32> de taille i  avec des valeurs entre -imax et + imax et controler les resultat
            let between = Uniform::<f32>::from(-fmax..fmax);
            let va: Vec<f32> = (0..i).map(|_| between.sample(&mut rng)).collect();
            let mut vb: Vec<f32> = (0..i).map(|_| between.sample(&mut rng)).collect();
            // reset half of vb to va
            vb[..(i / 2)].copy_from_slice(&va[..(i / 2)]);

            let easy_dist: u32 = va
                .iter()
                .zip(vb.iter())
                .map(|(a, b)| if a != b { 1 } else { 0 })
                .sum();
            let h_dist = DistHamming.eval(&va, &vb);
            let easy_dist = easy_dist as f32 / va.len() as f32;
            let j_exact = ((i / 2) as f32) / (i as f32);
            log::debug!(
                "test size {:?}  HammingDist {:.3e} easy : {:.3e} exact : {:.3e} ",
                i,
                h_dist,
                easy_dist,
                j_exact
            );
            if (easy_dist - h_dist).abs() > 1.0e-5 {
                println!(
                    " jhamming = {:?} , jexact = {:?}, j_easy : {:?}",
                    h_dist, j_exact, easy_dist
                );
                log::debug!("va = {:?}", va);
                log::debug!("vb = {:?}", vb);
                std::process::exit(1);
            }
            if (j_exact - h_dist).abs() > 1. / i as f32 + 1.0E-5 {
                println!(
                    " jhamming = {:?} , jexact = {:?}, j_easy : {:?}",
                    h_dist, j_exact, easy_dist
                );
                log::debug!("va = {:?}", va);
                log::debug!("vb = {:?}", vb);
                std::process::exit(1);
            }
        }
    } // end of test_hamming_f32

    #[cfg(feature = "stdsimd")]
    #[test]
    fn test_feature_simd() {
        init_log();
        log::info!("I have activated stdsimd");
    } // end of test_feature_simd

    #[test]
    #[cfg(feature = "simdeez_f")]
    fn test_feature_simdeez() {
        init_log();
        log::info!("I have activated simdeez");
    } // end of test_feature_simd

    #[test]
    #[cfg(feature = "half")]
    fn test_feature_half() {
        init_log();
        log::info!("I have activated half");
    } // end of test_feature_half

    #[test]
    #[cfg(feature = "half")]
    fn test_l1_f16() {
        use half::f16;

        let distl1 = DistL1;

        // Convert f32 values to f16
        let v1: Vec<f16> = vec![f16::from_f32(1.0), f16::from_f32(2.0), f16::from_f32(3.0)];
        let v2: Vec<f16> = vec![f16::from_f32(2.0), f16::from_f32(2.0), f16::from_f32(3.0)];

        let d1 = Distance::eval(&distl1, &v1, &v2);
        assert_eq!(d1, 1_f32);

        // Test with more values
        let v3: Vec<f16> = vec![
            f16::from_f32(1.5),
            f16::from_f32(-2.5),
            f16::from_f32(3.5),
            f16::from_f32(4.5),
        ];
        let v4: Vec<f16> = vec![
            f16::from_f32(1.0),
            f16::from_f32(-2.0),
            f16::from_f32(3.0),
            f16::from_f32(5.0),
        ];

        let d2 = distl1.eval(&v3, &v4);
        let expected = 0.5 + 0.5 + 0.5 + 0.5; // Sum of absolute differences
        assert!((d2 - expected).abs() < 1e-3);
    } // end of test_l1_f16

    #[test]
    #[cfg(feature = "half")]
    fn test_cosine_f16() {
        use half::f16;

        let distcos = DistCosine;

        // Simple orthogonal vectors test
        let v1: Vec<f16> = vec![f16::from_f32(1.0), f16::from_f32(0.0)];
        let v2: Vec<f16> = vec![f16::from_f32(0.0), f16::from_f32(1.0)];

        let d1 = Distance::eval(&distcos, &v1, &v2);
        assert!((d1 - 1.0).abs() < 1e-3); // Orthogonal vectors should have distance 1

        // Test with more realistic values
        let v3: Vec<f16> = vec![
            f16::from_f32(1.234),
            f16::from_f32(-1.678),
            f16::from_f32(1.367),
        ];
        let v4: Vec<f16> = vec![
            f16::from_f32(4.234),
            f16::from_f32(-6.678),
            f16::from_f32(10.367),
        ];

        let d2 = distcos.eval(&v3, &v4);

        // Calculate expected result manually
        let mut normv1 = 0.0;
        let mut normv2 = 0.0;
        let mut prod = 0.0;
        for i in 0..v3.len() {
            let a = v3[i].to_f32();
            let b = v4[i].to_f32();
            prod += a * b;
            normv1 += a * a;
            normv2 += b * b;
        }
        let expected = 1.0 - prod / (normv1 * normv2).sqrt();

        assert!((d2 - expected).abs() < 1e-3);
    } // end of test_cosine_f16

    #[test]
    #[cfg(feature = "half")]
    fn test_dot_f16() {
        use half::f16;

        let distdot = DistDot;

        // Create normalized vectors for testing
        // For DistDot, vectors should be L2 normalized
        let v1_raw = vec![1.0f32, 2.0f32, 3.0f32];
        let v2_raw = vec![4.0f32, 5.0f32, 6.0f32];

        // Normalize vectors
        let mut v1_norm = v1_raw.clone();
        let mut v2_norm = v2_raw.clone();
        l2_normalize(&mut v1_norm);
        l2_normalize(&mut v2_norm);

        // Convert to f16
        let v1: Vec<f16> = v1_norm.iter().map(|&x| f16::from_f32(x)).collect();
        let v2: Vec<f16> = v2_norm.iter().map(|&x| f16::from_f32(x)).collect();

        let d = distdot.eval(&v1, &v2);

        // Calculate expected result manually
        let dot_product: f32 = v1
            .iter()
            .zip(v2.iter())
            .map(|(a, b)| a.to_f32() * b.to_f32())
            .sum();
        let expected = 1.0 - dot_product;

        assert!((d - expected).abs() < 1e-3);
        assert!(expected >= 0.0); // Dot distance should be non-negative
    } // end of test_dot_f16

    #[test]
    #[cfg(feature = "half")]
    fn test_hellinger_f16() {
        use half::f16;

        // For Hellinger distance, vectors should be probability distributions
        // (non-negative and sum to 1)
        let length = 9;
        let mut p_data_f32 = Vec::with_capacity(length);
        let mut q_data_f32 = Vec::with_capacity(length);

        // Create probability distributions
        for _ in 0..length {
            p_data_f32.push(1.0 / length as f32);
            q_data_f32.push(1.0 / length as f32);
        }

        // Modify p slightly to make it different from q
        p_data_f32[0] -= 1.0 / (2 * length) as f32;
        p_data_f32[1] += 1.0 / (2 * length) as f32;

        // Convert to f16
        let p_data: Vec<f16> = p_data_f32.iter().map(|&x| f16::from_f32(x)).collect();
        let q_data: Vec<f16> = q_data_f32.iter().map(|&x| f16::from_f32(x)).collect();

        let dist = DistHellinger.eval(&p_data, &q_data);

        // Calculate expected result using the same formula as in test_hellinger
        let dist_exact_fn = |n: usize| -> f32 {
            let d1 = (4.0 - (6_f32).sqrt() - (2_f32).sqrt()) / n as f32;
            d1.sqrt() / (2_f32).sqrt()
        };
        let expected = dist_exact_fn(length);

        // Using a larger epsilon due to f16 precision limitations
        assert!((dist - expected).abs() < 1e-2);

        // Also verify that the result is in a reasonable range
        assert!(dist > 0.0);
        assert!(dist < 0.1); // The distance should be small for these similar distributions
    } // end of test_hellinger_f16

    #[test]
    #[cfg(feature = "half")]
    fn test_jeffreys_f16() {
        use half::f16;

        // For Jeffreys distance, vectors should be probability distributions
        let length = 9;
        let mut p_data_f32 = Vec::with_capacity(length);
        let mut q_data_f32 = Vec::with_capacity(length);

        // Create probability distributions
        for _ in 0..length {
            p_data_f32.push(1.0 / length as f32);
            q_data_f32.push(1.0 / length as f32);
        }

        // Modify distributions to make them different
        p_data_f32[0] -= 1.0 / (2 * length) as f32;
        p_data_f32[1] += 1.0 / (2 * length) as f32;
        q_data_f32[3] += 1.0 / (2 * length) as f32;

        // Convert to f16
        let p_data: Vec<f16> = p_data_f32.iter().map(|&x| f16::from_f32(x)).collect();
        let q_data: Vec<f16> = q_data_f32.iter().map(|&x| f16::from_f32(x)).collect();

        let dist_eval = DistJeffreys.eval(&p_data, &q_data);

        // Calculate expected result manually
        let mut dist_test = 0.0;
        for i in 0..length {
            let p = p_data_f32[i];
            let q = q_data_f32[i];
            dist_test += (p - q) * (p.max(M_MIN) / q.max(M_MIN)).ln();
        }

        assert!((dist_eval - dist_test).abs() < 1e-2); // Using a larger epsilon due to f16 precision
    } // end of test_jeffreys_f16

    #[test]
    #[cfg(feature = "half")]
    fn test_jensenshannon_f16() {
        use half::f16;

        // For Jensen-Shannon distance, vectors should be probability distributions
        let length = 9;
        let mut p_data_f32 = Vec::with_capacity(length);
        let mut q_data_f32 = Vec::with_capacity(length);

        // Create probability distributions
        for _ in 0..length {
            p_data_f32.push(1.0 / length as f32);
            q_data_f32.push(1.0 / length as f32);
        }

        // Modify distributions to make them different
        p_data_f32[0] -= 1.0 / (2 * length) as f32;
        p_data_f32[1] += 1.0 / (2 * length) as f32;
        q_data_f32[3] += 1.0 / (2 * length) as f32;

        // Convert to f16
        let p_data: Vec<f16> = p_data_f32.iter().map(|&x| f16::from_f32(x)).collect();
        let q_data: Vec<f16> = q_data_f32.iter().map(|&x| f16::from_f32(x)).collect();

        let dist_eval = DistJensenShannon.eval(&p_data, &q_data);

        // Calculate expected result manually
        let mut dist_test = 0.0;
        for i in 0..length {
            let p = p_data_f32[i];
            let q = q_data_f32[i];
            let mean_ab = 0.5 * (p + q);
            if p > 0.0 {
                dist_test += p * (p / mean_ab).ln();
            }
            if q > 0.0 {
                dist_test += q * (q / mean_ab).ln();
            }
        }
        let expected = (0.5 * dist_test).sqrt();

        assert!((dist_eval - expected).abs() < 1e-2); // Using a larger epsilon due to f16 precision
    } // end of test_jensenshannon_f16
} // end of module tests
