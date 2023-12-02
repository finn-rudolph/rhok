use core::cmp;
use std::{
    cmp::min,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use rand_xoshiro::{
    rand_core::{RngCore, SeedableRng},
    Xoshiro256PlusPlus,
};

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
pub fn pollard_rho(n: u64, k: u64, rng: &mut Xoshiro256PlusPlus) -> u64 {
    const BATCH_SIZE: u64 = 1 << 8;
    const LENGTH_LIMIT: u64 = 1 << 17;

    let mtg = Montgomery::new(n);

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

pub fn pollard_slow(n: u64, k: u64, rng: &mut Xoshiro256PlusPlus) -> u64 {
    const LENGTH_LIMIT: u64 = 1 << 18;

    let mtg = Montgomery::new(n);
    let _k = k << 1;

    loop {
        let mut x = rng.next_u64() % n;
        let mut y = x;

        for _ in 0..LENGTH_LIMIT {
            x = mtg.pow(x, _k) + 1;
            y = mtg.pow(mtg.pow(y, _k) + 1, _k) + 1;

            let mut d = x.wrapping_sub(y);
            if x < y {
                d = d.wrapping_add(n);
            }
            let g = gcd(d, n);
            if g != 1 && g != n {
                return g;
            }
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

// note: pollard_rho with BATCH_SIZE = 1 << 10 does not show improvement when
// M = 10 and some 2s and 3s are chosen (semiprimes 62).
// with 20/42 split, 2x2 and 2x3 outperforms 10x1.
// seems that when one prime factor is rather small, adding 2s and 3s helps
// something but WHY?

const K: [u64; 10] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

pub fn bench_single_rho() {
    const SAMPLES: usize = 10000;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos() as u64,
    );
    let mut avg = Duration::ZERO;

    for _ in 0..SAMPLES {
        let mut n = rng.next_u64() >> 2;
        while miller_rabin::miller_rabin(n) {
            n = rng.next_u64() >> 2;
        }
        let mut min_duration = Duration::from_secs(60);

        for k in K {
            let start = Instant::now();
            let factor = pollard_rho(n, k, &mut rng);
            min_duration = min(min_duration, start.elapsed());
        }

        avg += min_duration;
    }

    println!("{:?}", avg / SAMPLES as u32);
}
