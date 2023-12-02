use std::{ops::Range, time::Duration};

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rand_xoshiro::{
    rand_core::{RngCore, SeedableRng},
    Xoshiro256PlusPlus,
};

use rhok;

const K: Range<u64> = 1..221;
const BITS: [u32; 2] = [56, 62];

fn random_prime(bits: u32, rng: &mut Xoshiro256PlusPlus) -> u64 {
    let mut p = rng.next_u64() >> (64 - bits);
    while !rhok::miller_rabin::miller_rabin(p) {
        p = rng.next_u64() >> (64 - bits);
    }
    p
}

fn random_semiprime(bits: u32, rng: &mut Xoshiro256PlusPlus) -> u64 {
    random_prime(bits >> 1, rng) * random_prime(bits >> 1, rng)
}

fn bench_semiprime(c: &mut Criterion) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42u64 << 42);

    for bits in BITS {
        let mut group = c.benchmark_group(format!("single-semiprime-{}", bits));

        for k in K {
            group.bench_function(BenchmarkId::from_parameter(k), |b| {
                b.iter_batched(
                    || random_semiprime(bits, &mut rng),
                    |n| rhok::single::pollard_rho(black_box(n), black_box(k)),
                    criterion::BatchSize::SmallInput,
                )
            });
        }
    }
}

fn bench_general(c: &mut Criterion) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42u64 << 41);

    for bits in BITS {
        let mut group = c.benchmark_group(format!("single-general-{}", bits));

        for k in K {
            group.bench_function(BenchmarkId::from_parameter(k), |b| {
                b.iter_batched(
                    || {
                        random_prime(bits / 3, &mut rng)
                            * random_prime((2 * bits) / 3, &mut rng)
                    },
                    |n| rhok::single::pollard_rho(black_box(n), black_box(k)),
                    criterion::BatchSize::SmallInput,
                )
            });
        }
    }
}

criterion_group!(
    name = factorize_benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500)).sample_size(100);
    targets = bench_general, bench_semiprime
);

criterion_main!(factorize_benches);
