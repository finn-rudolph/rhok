use core::cmp;
use std::time::{Duration, Instant};

use rand_xoshiro::{
    rand_core::{RngCore, SeedableRng},
    Xoshiro256PlusPlus,
};
use rayon::prelude::*;

use crate::{miller_rabin, montgomery::Montgomery};

fn gcd(mut a: u64, mut b: u64) -> u64 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }

    let mut a_trz = a.trailing_zeros();
    let b_trz = b.trailing_zeros();
    let s = cmp::min(a_trz, b_trz);
    b >>= s;

    while a != 0 {
        a >>= a_trz;
        let d = b.wrapping_sub(a) as i64;
        a_trz = d.trailing_zeros();
        b = cmp::min(a, b);
        a = (if d < 0 { -d } else { d }) as u64;
    }

    b << s
}

// Pollard's Rho algorithm with Brent's optimization, taken from Sergey Slotin's
// book "Algorithms for Modern Hardware". (the gcd above is also from there)
fn pollard_rho(n: u64, k: u64) -> u64 {
    const BATCH_SIZE: u64 = 1 << 9;
    const LENGTH_LIMIT: u64 = 1 << 17;

    let mtg = Montgomery::new(n);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    loop {
        let mut x = rng.next_u64() % n;

        let mut l = BATCH_SIZE;
        while l <= LENGTH_LIMIT {
            let (y, mut q) = (x, 1);

            let mut i = 0;
            while i < l {
                for _ in 0..BATCH_SIZE {
                    x = mtg.pow(x, k << 1) + 1;
                    let mut d = x.wrapping_sub(y);
                    d = d.wrapping_add(if (d as i64) < 0 { n } else { 0 });
                    q = mtg.mul(q, d);
                }

                let g = gcd(q, n);
                if g != 1 && g != n {
                    return g;
                }

                i += BATCH_SIZE;
            }

            l <<= 1;
        }
    }
}

fn random_prime(bits: u32, rng: &mut Xoshiro256PlusPlus) -> u64 {
    let mut p = rng.next_u64() >> (64 - bits);
    while !miller_rabin::miller_rabin(p) {
        p = rng.next_u64() >> (64 - bits);
    }
    p
}

pub fn bench_single_rho() {
    const SAMPLES: usize = 1000;

    print!("[");
    (2..10)
        .map(|k| {
            let mut duration_sum = (0..SAMPLES)
                .into_par_iter()
                .fold(
                    || Duration::ZERO,
                    |sum, i| {
                        let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                            (k << 42) ^ i as u64,
                        );
                        let n = random_prime(31, &mut rng)
                            * random_prime(31, &mut rng);
                        let start = Instant::now();
                        let factor = pollard_rho(n, k);
                        sum + start.elapsed()
                    },
                )
                .sum::<Duration>();
            (k, duration_sum / SAMPLES as u32)
        })
        .collect::<Vec<(u64, Duration)>>()
        .iter()
        .for_each(|(k, t)| print!("{},", t.as_nanos()));
    print!("]");
}
