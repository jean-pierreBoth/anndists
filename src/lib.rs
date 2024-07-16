#![cfg_attr(feature = "stdsimd", feature(portable_simd))]
//
// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::Builder;

use lazy_static::lazy_static;

pub mod dist;

pub mod prelude;

lazy_static! {
    static ref LOG: u64 = init_log();
}

// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    1
}
