use core::cmp;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use rand::{seq::SliceRandom, thread_rng, Rng};
use rand_xoshiro::{
    rand_core::{RngCore, SeedableRng},
    Xoshiro256PlusPlus,
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{
    miller_rabin::{self, miller_rabin},
    montgomery::Montgomery,
};

pub fn gcd(mut a: u64, mut b: u64) -> u64 {
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
    const BATCH_SIZE: u64 = 1 << 6;
    const LENGTH_LIMIT: u64 = 1 << 17;

    let mtg = Montgomery::new(n);
    let one = mtg.to_montgomery_space(1);
    let _k = k << 1;

    loop {
        let mut x = rng.next_u64() % n;

        let mut l = BATCH_SIZE;
        while l <= LENGTH_LIMIT {
            let (y, mut q) = (x, 1);

            let mut i = 0;
            while i < l {
                for _ in 0..BATCH_SIZE {
                    x = mtg.add(mtg.pow(x, _k), one);
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

pub fn pollard_rho_iteration_count(
    n: u64,
    k: u64,
    rng: &mut Xoshiro256PlusPlus,
) -> usize {
    const BATCH_SIZE: u64 = 1 << 6;
    const LENGTH_LIMIT: u64 = 1 << 17;

    let mtg = Montgomery::new(n);
    let one = mtg.to_montgomery_space(1);
    let _k = k << 1;
    let mut iterations: usize = 0;

    loop {
        let mut x = rng.next_u64() % n;

        let mut l = BATCH_SIZE;
        while l <= LENGTH_LIMIT {
            let (y, mut q) = (x, 1);

            let mut i = 0;
            while i < l {
                for _ in 0..BATCH_SIZE {
                    iterations += 1;
                    x = mtg.add(mtg.pow(x, _k), one);
                    let mut d = x.wrapping_sub(y);
                    d = d.wrapping_add(if (d as i64) < 0 { n } else { 0 });
                    q = mtg.mul(q, d);
                }

                let g = gcd(q, n);
                if g != 1 {
                    return iterations;
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
    let one = mtg.to_montgomery_space(1);
    let _k = k << 1;

    loop {
        let mut x = rng.next_u64() % n;
        let mut y = x;

        for _ in 0..LENGTH_LIMIT {
            x = mtg.add(mtg.pow(x, _k), one);
            y = mtg.add(mtg.pow(mtg.add(mtg.pow(y, _k), one), _k), one);

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

pub fn fn_iteration_count(
    v: &Vec<usize>,
    rng: &mut Xoshiro256PlusPlus,
) -> usize {
    const LENGTH_LIMIT: u64 = 1 << 18;

    let mut iterations: usize = 0;

    loop {
        let mut x = rng.next_u64() as usize % v.len();
        let mut y = x;

        for _ in 0..LENGTH_LIMIT {
            iterations += 1;
            x = v[x];
            y = v[v[y]];

            if y == x {
                return iterations;
            }
        }

        assert!(false);
    }
}

pub fn pollard_slow_iteration_count(
    p: u64,
    k: u64,
    rng: &mut Xoshiro256PlusPlus,
) -> usize {
    let mtg = Montgomery::new(p);
    let mut iterations: usize = 0;

    let mut x = rng.next_u64() % p;
    let mut y = x;
    let _k = k << 1;

    loop {
        iterations += 1;
        x = mtg.strict(mtg.add(mtg.pow(x, _k), mtg.one()));
        y = mtg.strict(
            mtg.add(mtg.pow(mtg.add(mtg.pow(y, _k), mtg.one()), _k), mtg.one()),
        );
        if x == y {
            return iterations;
        }
    }
}

fn _random_prime(bits: u32, rng: &mut Xoshiro256PlusPlus) -> u64 {
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

const K: u64 = 3;

pub fn bench_single_rho() {
    const A: usize = 1 << 42;
    const B: usize = (1 << 42) + (1 << 14);
    const ITER: usize = 100;

    let samples: Vec<f64> = (A..=B)
        .into_par_iter()
        .map(|n| {
            if miller_rabin(n as u64) {
                return 0.0;
            }

            let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos() as u64,
            );

            let start = Instant::now();
            for _ in 0..ITER {
                pollard_rho(n as u64, K, &mut rng);
            }

            start.elapsed().as_nanos() as f64
                / (ITER as f64 * (n as f64).sqrt())
        })
        .filter(|x| *x != 0.0)
        .collect();

    println!(
        "{:?}",
        (samples.iter().fold(0.0, |acc, elem| acc + elem)
            / samples.len() as f64)
    );
}
