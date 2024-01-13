use core::cmp;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

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
    const BATCH_SIZE: u64 = 1 << 10;
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
                if g != 1 {
                    return g;
                }

                i += BATCH_SIZE;
            }

            l <<= 1;
        }
    }
}

pub fn bench_single_rho() {
    const A: usize = 1 << 61;
    const B: usize = (1 << 61) + (1 << 14);
    const ITERATIONS: usize = 10;

    for k1 in 1..17 {
        for k2 in 1..17 {
            let mut total_time = Duration::ZERO;
            let mut num_composites = 0;

            for n in A..=B {
                if miller_rabin(n as u64) {
                    continue;
                }

                num_composites += 1;

                let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_nanos() as u64,
                );

                for _ in 0..ITERATIONS {
                    let mut min_time;

                    let mut start = Instant::now();
                    pollard_rho(n as u64, k1, &mut rng);
                    min_time = start.elapsed();

                    start = Instant::now();
                    pollard_rho(n as u64, k2, &mut rng);
                    min_time = min_time.min(start.elapsed());

                    total_time += min_time;
                }
            }

            println!(
                "k1 = {:<5} | k2 = {:<5} | {:?}",
                k1,
                k2,
                total_time / (num_composites * ITERATIONS) as u32
            );
        }
    }
}
