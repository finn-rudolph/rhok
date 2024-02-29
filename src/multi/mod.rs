mod pollard_rho;

use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use rug::{integer::IsPrime, rand::RandState, Integer};

use self::pollard_rho::pollard_rho;

const MIN_BITS: u32 = 22;
const MAX_BITS: u32 = 128;
const TOTAL_BITS: u32 = 192;

fn prime_with_bits(
    min_bits: u32,
    max_bits: u32,
    rng: &mut RandState,
) -> Integer {
    let mut p: Integer = Integer::random_bits(max_bits, rng).into();
    while p.significant_bits() < min_bits
        || p.is_probably_prime(32) == IsPrime::No
    {
        p = Integer::random_bits(max_bits, rng).into();
    }
    p
}

pub fn gen_test_numbers(samples: usize) -> Vec<Integer> {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos(),
    ));
    let mut test_numbers: Vec<Integer> = Vec::new();
    while test_numbers.len() < samples {
        let mut n = prime_with_bits(MIN_BITS, MIN_BITS, &mut rng);
        while n.significant_bits() < TOTAL_BITS {
            n *= prime_with_bits(MIN_BITS, MIN_BITS, &mut rng);
        }
        test_numbers.push(n);
    }
    test_numbers
}

pub fn measure(k: &Vec<u64>, test_numbers: &Vec<Integer>) -> (f64, usize) {
    let mut sum = Duration::ZERO;
    let mut outliers: usize = 0;
    let mut rng = RandState::new();

    for n in test_numbers {
        let mut min_time = Duration::from_secs(42);

        for k_i in k {
            let start = Instant::now();
            match pollard_rho(n, *k_i, &mut rng) {
                Some(_) => min_time = min_time.min(start.elapsed()),
                None => {
                    outliers += 1;
                    min_time = Duration::ZERO;
                    break;
                }
            }
        }

        sum += min_time;
    }

    (
        sum.as_nanos() as f64 / (test_numbers.len() - outliers) as f64,
        outliers,
    )
}
