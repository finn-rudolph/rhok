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

const SAMPLES: usize = 1 << 16;

#[derive(Clone, Copy, PartialEq, Eq)]
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
    k_min: u64,
    k_max: u64,
    source: Source,
    raw: bool,
    k: &mut Vec<u64>,
    j: usize,
) {
    if j == k.len() {
        if !raw {
            for k_j in k.iter() {
                print!("{:<5}", k_j);
            }
        }

        let val = match source {
            Source::Real => {
                let mut samples: Vec<Duration> = (0..SAMPLES)
                    .into_par_iter()
                    .map(|_| {
                        let mut rng = thread_rng();
                        let n = random_prime(21, &mut rng)
                            * random_prime(41, &mut rng);

                        let mut min_time = Duration::from_secs(42);
                        for k_j in k.iter() {
                            let start = Instant::now();
                            pollard_rho::pollard_rho(n, *k_j, &mut rng);
                            let time_needed = start.elapsed();
                            if time_needed >= Duration::from_millis(200) {
                                return Duration::ZERO;
                            }
                            min_time = min_time.min(time_needed);
                        }

                        min_time
                    })
                    .collect();

                samples.sort();
                let mut start = 0;
                while samples[start] == Duration::ZERO {
                    start += 1;
                }

                print!(
                    "{:<20}{:<20}",
                    (start as f64 / SAMPLES as f64) * 100.0,
                    start
                );

                samples[start..].iter().sum::<Duration>().as_nanos() as f64
                    / (SAMPLES - start) as f64
            }

            Source::Formula => formula::expected_time(k_min, k_max, k, 0, 0.0),
        };
        println!("{}", val);

        return;
    }

    k[j] = if j > 0 { k[j - 1] } else { k_min };
    while k[j] <= k_max {
        iterate_k_cartesian_product(k_min, k_max, source, raw, k, j + 1);
        k[j] += 1;
    }
}

// CLI usage: [--formula | --real] [#machines] [minimal k] [maximal k]
// If `--raw` is written after this, only the values without additional
// information are printed.

fn main() {
    // Only use half of the available threads for better measurement accuracy.
    rayon::ThreadPoolBuilder::new()
        .num_threads(10)
        .build_global()
        .unwrap();

    let args: Vec<String> = env::args().collect();
    assert!(5 <= args.len() && args.len() <= 6);

    let source = match args[1].as_str() {
        "--formula" => Source::Formula,
        "--real" => Source::Real,
        _ => unreachable!(),
    };

    let machines: usize = args[2].parse().unwrap();
    let k_min: u64 = args[3].parse().unwrap();
    let k_max: u64 = args[4].parse().unwrap();

    let raw = if args.len() == 6 {
        assert!(args[5].as_str() == "--raw");
        true
    } else {
        for i in 1..=machines {
            print!("{:<5}", format!("k{}", i));
        }

        if source == Source::Real {
            print!("{:<20}{:<20}", "% outliers", "# outliers");
        }

        println!("value");

        false
    };

    assert!(machines > 0);
    assert!(k_min > 0);
    assert!(k_min <= k_max);

    let mut k = vec![0u64; machines];

    iterate_k_cartesian_product(k_min, k_max, source, raw, &mut k, 0);
}
