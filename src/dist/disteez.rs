//! simdeez distance implementations
//!
//!

#![cfg(feature = "simdeez_f")]
use simdeez::*;
use simdeez::{prelude::*, simd_runtime_generate};

use super::distances::M_MIN;

simd_runtime_generate!(
    pub(super) fn distance_l1_f32_simdeez(va: &[f32], vb: &[f32]) -> f32 {
        assert_eq!(va.len(), vb.len());
        //
        let mut dist: f32;
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        let mut dist_simd = S::Vf32::zeroes();
        //
        while a.len() >= S::Vf32::WIDTH {
            let xa = S::Vf32::load_from_slice(a);
            let xb = S::Vf32::load_from_slice(b);
            let delta = S::Vf32::abs(xa - xb);
            dist_simd += delta;
            //
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        // horizontal add
        dist = S::Vf32::horizontal_add(dist_simd);
        // remaining
        for i in 0..a.len() {
            //        log::debug!("distance_l1_f32, i {:?} len {:?} ", i, a.len());
            dist += (a[i] - b[i]).abs();
        }
        assert!(dist >= 0.);
        dist
    }
);

//========================================================================

simd_runtime_generate!(
    pub(super) fn distance_l2_f32_simdeez(va: &[f32], vb: &[f32]) -> f32 {
        //
        assert_eq!(va.len(), vb.len());
        //
        let mut dist: f32;
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        let mut dist_simd = S::Vf32::zeroes();

        while a.len() >= S::Vf32::WIDTH {
            let xa = S::Vf32::load_from_slice(a);
            let xb = S::Vf32::load_from_slice(b);
            let mut delta = xa - xb;
            delta *= delta;
            dist_simd += delta;
            // shift
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        //
        dist = S::Vf32::horizontal_add(dist_simd);
        // remaining
        for i in 0..a.len() {
            dist += (a[i] - b[i]) * (a[i] - b[i]);
        }
        assert!(dist >= 0.);
        dist.sqrt()
    }
);

//======================================================================

simd_runtime_generate!(
    pub(super) fn distance_dot_f32_simdeez(va: &[f32], vb: &[f32]) -> f32 {
        //
        assert_eq!(va.len(), vb.len());
        //
        let mut dot: f32;
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        let mut dot_simd = S::Vf32::zeroes();
        //
        while a.len() >= S::Vf32::WIDTH {
            let xa = S::Vf32::load_from_slice(a);
            let xb = S::Vf32::load_from_slice(b);
            let delta = xa * xb;
            dot_simd += delta;
            // shift
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        dot = S::Vf32::horizontal_add(dot_simd);
        //
        for i in 0..a.len() {
            dot += a[i] * b[i];
        }
        assert!(dot <= 1.000002);
        (1. - dot).max(0.)
    }
);

//============================================================================

simd_runtime_generate!(
    pub(super) fn distance_hellinger_f32_simdeez(va: &[f32], vb: &[f32]) -> f32 {
        assert_eq!(va.len(), vb.len());
        //
        log::debug!("in disteez::distance_hellinger_f32_simdeez");
        //
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        let mut dist_simd = S::Vf32::zeroes();
        let mut dist: f32;
        //
        while a.len() >= S::Vf32::WIDTH {
            let xa = S::Vf32::load_from_slice(a);
            let xb = S::Vf32::load_from_slice(b);
            let prod = xa * xb;
            let prod_s = S::Vf32::sqrt(prod);
            dist_simd += prod_s;
            // shift
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        dist = S::Vf32::horizontal_add(dist_simd);

        for i in 0..a.len() {
            dist += a[i].sqrt() * b[i].sqrt();
        }
        assert!(1. - dist >= -0.000001);
        dist = (1. - dist).max(0.).sqrt();
        dist
    }
);

//=============================================================================

simd_runtime_generate!(
    pub(super) fn distance_jeffreys_f32_simdeez(va: &[f32], vb: &[f32]) -> f32 {
        //
        log::debug!("in disteez::distance_jeffreys_f32_simdeez");
        //
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        let mut logslice = Vec::<f32>::with_capacity(S::Vf32::WIDTH);
        let mut dist_simd = S::Vf32::zeroes();
        let mut dist: f32;
        //
        while a.len() >= S::Vf32::WIDTH {
            let xa = S::Vf32::load_from_slice(a);
            let xb = S::Vf32::load_from_slice(b);
            let delta = xa - xb;
            for j in 0..S::Vf32::WIDTH {
                // take care of zeros!
                logslice.push((xa[j].max(M_MIN) / xb[j].max(M_MIN)).ln());
            }
            let prod_s = delta * S::Vf32::load_from_slice_exact(logslice.as_slice()).unwrap();
            dist_simd += prod_s;
            logslice.clear();
            // shift
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        dist = S::Vf32::horizontal_add(dist_simd);
        // remaining
        for i in 0..a.len() {
            if vb[i] > 0. {
                dist += (a[i] - b[i]) * (a[i].max(M_MIN) / b[i].max(M_MIN)).ln();
            }
        }
        dist
    }
);

simd_runtime_generate!(
    pub(super) fn distance_hamming_i32_simdeez(va: &[i32], vb: &[i32]) -> f32 {
        assert_eq!(va.len(), vb.len());
        //
        log::debug!(" in disteez::distance_hamming_i32_simdeez");
        //
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        //
        let mut dist: i32;
        let mut simd_res: Vec<i32> = (0..S::Vi32::WIDTH).map(|_| 0).collect();
        //
        let mut dist_simd = S::Vi32::zeroes();
        while a.len() >= S::Vf32::WIDTH {
            let xa = S::Vi32::load_from_slice(a);
            let xb = S::Vi32::load_from_slice(b);
            let delta = S::Vi32::cmp_neq(xa, xb);
            dist_simd += delta;
            // shift
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        // get the sum of value in dist
        dist_simd.copy_to_slice(&mut simd_res);

        dist = simd_res.into_iter().sum();
        //
        // Beccause simd returns 0xFFFF... when neq true and 0 else
        dist = -dist;
        // add the residue
        for i in 0..a.len() {
            dist += if a[i] != b[i] { 1 } else { 0 }
        }
        dist as f32 / va.len() as f32
    }
);

macro_rules! unsigned_hamming_simdeez {
    ($name:ident, $unsigned:ty, $signed:ty, $vector:ident) => {
        simd_runtime_generate!(
            #[allow(dead_code)] // Also generates scalar and generic entry points.
            pub(super) fn $name(va: &[$unsigned], vb: &[$unsigned]) -> f32 {
                assert_eq!(va.len(), vb.len());

                // Each lane gains at most one per vector; bound even the i16 counters.
                let block_size = i16::MAX as usize * S::$vector::WIDTH;
                let mut dist = 0usize;
                for (a, b) in va.chunks(block_size).zip(vb.chunks(block_size)) {
                    let mut a_chunks = a.chunks_exact(S::$vector::WIDTH);
                    let mut b_chunks = b.chunks_exact(S::$vector::WIDTH);
                    let mut counts = S::$vector::zeroes();
                    for (a_chunk, b_chunk) in a_chunks.by_ref().zip(b_chunks.by_ref()) {
                        // SAFETY: each chunk has WIDTH elements. Same-width signed and
                        // unsigned integers share layout and accept every bit pattern;
                        // the loads do not require SIMD alignment.
                        let (a_vec, b_vec) = unsafe {
                            (
                                S::$vector::load_from_ptr_unaligned(
                                    a_chunk.as_ptr().cast::<$signed>(),
                                ),
                                S::$vector::load_from_ptr_unaligned(
                                    b_chunk.as_ptr().cast::<$signed>(),
                                ),
                            )
                        };
                        // A comparison mask is -1 for a mismatch and 0 for a match.
                        counts -= a_vec.cmp_neq(b_vec);
                    }
                    dist += counts.horizontal_add() as usize;
                    dist += a_chunks
                        .remainder()
                        .iter()
                        .zip(b_chunks.remainder())
                        .filter(|(a, b)| a != b)
                        .count();
                }
                dist as f32 / va.len() as f32
            }
        );
    };
}

unsigned_hamming_simdeez!(distance_hamming_u16_simdeez, u16, i16, Vi16);
unsigned_hamming_simdeez!(distance_hamming_u32_simdeez, u32, i32, Vi32);
unsigned_hamming_simdeez!(distance_hamming_u64_simdeez, u64, i64, Vi64);

// special implementation for f64 exclusively in the context of SuperMinHash algorithm
simd_runtime_generate!(
    pub(super) fn distance_hamming_f64(va: &[f64], vb: &[f64]) -> f32 {
        assert_eq!(va.len(), vb.len());
        //
        log::debug!(" in disteez::distance_hamming_f64");
        //
        let mut a = &va[..va.len()];
        let mut b = &vb[..va.len()];
        let mut simd_res: Vec<i64> = (0..S::Vf64::WIDTH).map(|_| 0).collect();
        //
        let mut dist_simd = S::Vi64::zeroes();
        //    log::debug!("initial simd_res : {:?}", dist_simd);
        while a.len() >= S::Vf64::WIDTH {
            let xa = S::Vf64::load_from_slice(a);
            let xb = S::Vf64::load_from_slice(b);
            let delta = S::Vf64::cmp_neq(xa, xb);
            let delta_i = delta.bitcast_i64();
            //        log::debug!("delta_i : , {:?}", delta_i);
            // cast to i64 to transform the 0xFFFFF.... to -1
            dist_simd += delta_i;
            // shift
            a = &a[S::Vf64::WIDTH..];
            b = &b[S::Vf64::WIDTH..];
        }
        // get the sum of value in dist
        //    log::trace!("simd_res : {:?}", dist_simd);
        dist_simd.copy_to_slice(&mut simd_res);

        // cmp_neq returns 0xFFFFFFFFFF if true and 0 else, we need to transform 0xFFFFFFF... to 1
        simd_res.iter_mut().for_each(|x| *x = -*x);
        //    log::debug!("simd_res : {:?}", simd_res);
        let mut dist: i64 = simd_res.into_iter().sum();
        // Beccause simd returns 0xFFFF... when neq true and 0 else
        // add the residue
        for i in 0..a.len() {
            dist += if a[i] != b[i] { 1 } else { 0 }
        }
        (dist as f64 / va.len() as f64) as f32
    }
);

//=======================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::*;

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    use rand::distr::{Distribution, Uniform};

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    fn init_log() -> u64 {
        let mut builder = env_logger::Builder::from_default_env();
        let _ = builder.is_test(true).try_init();
        println!("\n ************** initializing logger *****************\n");
        1
    }

    fn check_unsigned_hamming_backend<S: Simd>() {
        macro_rules! check_width {
            ($unsigned:ty, $vector:ident, $kernel:ident) => {
                let check = |a: &[$unsigned], b: &[$unsigned]| {
                    let count = a.iter().zip(b).filter(|(a, b)| a != b).count();
                    let want = count as f32 / a.len() as f32;
                    let got = $kernel::<S>(a, b);
                    let public = DistHamming.eval(a, b);
                    if want.is_nan() {
                        assert!(got.is_nan());
                        assert!(public.is_nan());
                    } else {
                        assert_eq!(got.to_bits(), want.to_bits());
                        assert_eq!(got.to_bits(), public.to_bits());
                    }
                };
                let block_boundary = i16::MAX as usize * S::$vector::WIDTH;
                for len in (0..=65).chain([
                    block_boundary - 1,
                    block_boundary,
                    block_boundary + 1,
                    block_boundary + S::$vector::WIDTH + 1,
                ]) {
                    let a = vec![<$unsigned>::MAX; len + 1];
                    let mut b = a.clone();
                    check(&a[1..], &b[1..]);
                    for (i, value) in b[1..].iter_mut().enumerate() {
                        *value ^= (1 as $unsigned) << (i % <$unsigned>::BITS as usize);
                    }
                    check(&a[1..], &b[1..]);
                    for i in (1..b.len()).step_by(3) {
                        b[i] = a[i];
                    }
                    check(&a[1..], &b[1..]);
                }
            };
        }

        check_width!(u16, Vi16, distance_hamming_u16_simdeez_generic);
        check_width!(u32, Vi32, distance_hamming_u32_simdeez_generic);
        check_width!(u64, Vi64, distance_hamming_u64_simdeez_generic);
    }

    #[test]
    fn test_unsigned_hamming_scalar() {
        check_unsigned_hamming_backend::<engines::scalar::Scalar>();
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_unsigned_hamming_neon() {
        if std::arch::is_aarch64_feature_detected!("neon") {
            check_unsigned_hamming_backend::<engines::neon::Neon>();
        }
    }

    #[test]
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    fn test_unsigned_hamming_sse2() {
        if is_x86_feature_detected!("sse2") {
            check_unsigned_hamming_backend::<engines::sse2::Sse2>();
        }
    }

    #[test]
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    fn test_unsigned_hamming_sse41() {
        if is_x86_feature_detected!("sse4.1") {
            check_unsigned_hamming_backend::<engines::sse41::Sse41>();
        }
    }

    #[test]
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    fn test_unsigned_hamming_avx2() {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            check_unsigned_hamming_backend::<engines::avx2::Avx2>();
        }
    }

    #[test]
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    fn test_unsigned_hamming_avx512() {
        if is_x86_feature_detected!("avx512f")
            && is_x86_feature_detected!("avx512bw")
            && is_x86_feature_detected!("avx512dq")
        {
            check_unsigned_hamming_backend::<engines::avx512::Avx512>();
        }
    }

    #[test]
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    fn test_unsigned_hamming_wasm() {
        check_unsigned_hamming_backend::<engines::wasm32::Wasm>();
    }

    #[test]
    fn test_simdeez_hamming_i32() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            init_log();
            log::info!("running test_simd_hamming_i32 for avx2");
            //
            let size_test = 500;
            let imax = 3;
            let mut rng = rand::rng();
            for i in 4..size_test {
                // generer 2 va et vb s des vecteurs<i32> de taille i  avec des valeurs entre -imax et + imax et controler les resultat
                let between = Uniform::<i32>::try_from(-imax..imax).unwrap();
                let va: Vec<i32> = (0..i)
                    .into_iter()
                    .map(|_| between.sample(&mut rng))
                    .collect();
                let vb: Vec<i32> = (0..i)
                    .into_iter()
                    .map(|_| between.sample(&mut rng))
                    .collect();
                let simd_dist = distance_hamming_i32_simdeez(&va, &vb);

                let easy_dist: u32 = va
                    .iter()
                    .zip(vb.iter())
                    .map(|(a, b)| if a != b { 1 } else { 0 })
                    .sum();
                let easy_dist = easy_dist as f32 / va.len() as f32;
                log::debug!(
                    "test size {:?} simd  exact = {:?} {:?}",
                    i,
                    simd_dist,
                    easy_dist
                );
                if (easy_dist - simd_dist).abs() > 1.0e-5 {
                    println!(" jsimd = {:?} , easy dist = {:?}", simd_dist, easy_dist);
                    println!("va = {:?}", va);
                    println!("vb = {:?}", vb);
                    std::process::exit(1);
                }
            }
        } // cfg
    } // end of test_simdeez_hamming_i32

    #[test]
    fn test_simdeez_hamming_f64() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            init_log();
            log::info!("running test_simd_hamming_f64 for avx2");
            //
            let size_test = 500;
            let fmax: f64 = 3.;
            let mut rng = rand::rng();
            for i in 300..size_test {
                // generer 2 va et vb s des vecteurs<i32> de taille i  avec des valeurs entre -imax et + imax et controler les resultat
                let between = Uniform::<f64>::try_from(-fmax..fmax).unwrap();
                let va: Vec<f64> = (0..i)
                    .into_iter()
                    .map(|_| between.sample(&mut rng))
                    .collect();
                let mut vb: Vec<f64> = (0..i)
                    .into_iter()
                    .map(|_| between.sample(&mut rng))
                    .collect();
                // reset half of vb to va
                vb[..(i / 2)].copy_from_slice(&va[..(i / 2)]);
                let simd_dist = distance_hamming_f64(&va, &vb);

                let j_exact = ((i / 2) as f32) / (i as f32);
                let easy_dist: u32 = va
                    .iter()
                    .zip(vb.iter())
                    .map(|(a, b)| if a != b { 1 } else { 0 })
                    .sum();
                let h_dist = DistHamming.eval(&va, &vb);
                let easy_dist = easy_dist as f32 / va.len() as f32;
                log::debug!(
                    "test size {:?} simd  = {:.3e} HammingDist {:.3e} easy : {:.3e} exact : {:.3e} ",
                    i,
                    simd_dist,
                    h_dist,
                    easy_dist,
                    0.5
                );
                if (easy_dist - simd_dist).abs() > 1.0e-5 {
                    println!(" jsimd = {:?} , jexact = {:?}", simd_dist, easy_dist);
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
        } // cfg
    } // end of test_simd_hamming_f64
} // end of mod tests
