use std::time::Duration;

mod miller_rabin;
mod montgomery;
mod pollard_rho;

pub use miller_rabin::miller_rabin;

pub fn measure(k: &Vec<u64>) -> f64 {
    todo!()
}
