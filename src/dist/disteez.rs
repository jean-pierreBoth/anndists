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
            let xa = S::Vf32::load_from_slice(&a);
            let xb = S::Vf32::load_from_slice(&b);
            let delta = S::Vf32::abs(xa - xb);
            dist_simd += delta;
            //
            a = &a[S::Vf32::WIDTH..];
            b = &b[S::Vf32::WIDTH..];
        }
        // horizontal add
        dist = S::Vf32::horizontal_add(dist_simd);
        // remaining
        for i in 0..va.len() {
            //        log::debug!("distance_l1_f32, i {:?} len {:?} ", i, va.len());
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
            let xa = S::Vf32::load_from_slice(&a);
            let xb = S::Vf32::load_from_slice(&b);
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
        for i in 0..va.len() {
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
            let xa = S::Vf32::load_from_slice(&a);
            let xb = S::Vf32::load_from_slice(&b);
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
            let xa = S::Vf32::load_from_slice(&a);
            let xb = S::Vf32::load_from_slice(&b);
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
            let xa = S::Vf32::load_from_slice(&a);
            let xb = S::Vf32::load_from_slice(&b);
            let delta = xa - xb;
            for j in 0..S::Vf32::WIDTH {
                // take care of zeros!
                logslice.push((xa[j].max(M_MIN) / xb[j].max(M_MIN)).ln());
            }
            let prod_s = delta * S::Vf32::load_from_slice_exact(&logslice.as_slice()).unwrap();
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
            let xa = S::Vi32::load_from_slice(&a);
            let xb = S::Vi32::load_from_slice(&b);
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
            let xa = S::Vf64::load_from_slice(&a);
            let xb = S::Vf64::load_from_slice(&b);
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

    use rand::distr::{Distribution, Uniform};

    fn init_log() -> u64 {
        let mut builder = env_logger::Builder::from_default_env();
        let _ = builder.is_test(true).try_init();
        println!("\n ************** initializing logger *****************\n");
        return 1;
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
                let simd_dist = distance_hamming_i32_simdeez(&va, &vb) as f32;

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
                for i in 0..i / 2 {
                    vb[i] = va[i];
                }
                let simd_dist = distance_hamming_f64(&va, &vb) as f32;

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
