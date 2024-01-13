mod formula;
mod miller_rabin;
mod montgomery;
mod pollard_rho;

use std::{
    env,
    time::{Duration, Instant},
};

use rand::{thread_rng, RngCore};
use rayon::prelude::*;

use crate::miller_rabin::miller_rabin;

const SAMPLES: usize = 1 << 20;

#[derive(Clone, Copy)]
enum Source {
    Formula,
    Real,
}

fn random_prime(bits: u32, rng: &mut dyn RngCore) -> u64 {
    loop {
        let p = rng.next_u64() >> (64 - bits);
        if p != 2 && miller_rabin(p) {
            return p;
        }
    }
}

fn iterate_k_cartesian_product(
    machines: usize,
    k_min: u64,
    k_max: u64,
    source: Source,
    k: &mut Vec<u64>,
    j: usize,
) {
    if j == machines {
        for k_j in k.iter() {
            print!("{:<5}", k_j);
        }

        println!(
            " | {}",
            match source {
                Source::Real => {
                    let total_time: Duration = (0..SAMPLES)
                        .into_par_iter()
                        .map(|_| {
                            let mut rng = thread_rng();
                            let n = random_prime(27, &mut rng)
                                * random_prime(31, &mut rng);

                            let mut min_time = Duration::from_secs(42);
                            for k_j in k.iter() {
                                let start_time = Instant::now();
                                pollard_rho::pollard_rho(n, *k_j, &mut rng);
                                min_time = min_time.min(start_time.elapsed());
                            }

                            min_time
                        })
                        .filter(|x| *x != Duration::ZERO)
                        .sum();

                    total_time.as_nanos() as f64 / SAMPLES as f64
                }
                Source::Formula =>
                    formula::expected_time(machines, k_min, k_max, k, 0, 0.0),
            }
        );

        return;
    }

    k[j] = if j > 0 { k[j - 1] } else { k_min };
    while k[j] <= k_max {
        iterate_k_cartesian_product(machines, k_min, k_max, source, k, j + 1);
        k[j] += 1;
    }
}

fn main() {
    // Only use half of the available threads for better measurement accuracy.
    rayon::ThreadPoolBuilder::new()
        .num_threads(10)
        .build_global()
        .unwrap();

    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 5);

    let source = match args[1].as_str() {
        "--formula" => Source::Formula,
        "--real" => Source::Real,
        _ => unreachable!(),
    };

    let machines: usize = args[2].parse().unwrap();
    let k_min: u64 = args[3].parse().unwrap();
    let k_max: u64 = args[4].parse().unwrap();

    assert!(machines > 0);
    assert!(k_min > 0);
    assert!(k_min <= k_max);

    let mut k = vec![0u64; machines];

    iterate_k_cartesian_product(machines, k_min, k_max, source, &mut k, 0);
}
