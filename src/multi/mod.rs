mod pollard_rho;

use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use once_cell::sync::Lazy;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rug::{integer::IsPrime, rand::RandState, Integer};

use self::pollard_rho::pollard_rho;

const SAMPLES: usize = 1 << 10;
const MIN_BITS: u32 = 22;
const MAX_BITS: u32 = 128;
const TOTAL_BITS: u32 = 192;

static TEST_NUMBERS: Lazy<Vec<Integer>> =
    Lazy::new(|| gen_test_numbers(SAMPLES));

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

fn gen_test_numbers(samples: usize) -> Vec<Integer> {
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
            n *= prime_with_bits(
                MIN_BITS,
                (Integer::from(MAX_BITS - MIN_BITS).random_below(&mut rng)
                    + MIN_BITS)
                    .to_u32()
                    .unwrap(),
                &mut rng,
            );
        }
        test_numbers.push(n);
    }

    test_numbers
}

pub fn measure(k: &Vec<u64>) -> f64 {
    let sum: Duration = TEST_NUMBERS
        .par_iter()
        .fold(
            || Duration::ZERO,
            |sum, n| {
                let mut min_time = Duration::from_secs(42);
                let mut rng = RandState::new();

                for k_i in k {
                    let start = Instant::now();
                    pollard_rho(n, *k_i, &mut rng);
                    min_time = min_time.min(start.elapsed());
                }

                sum + min_time
            },
        )
        .sum();

    sum.as_nanos() as f64 / TEST_NUMBERS.len() as f64
}
