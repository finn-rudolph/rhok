use std::{ops::Range, time::Duration};

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rug::{integer::IsPrime, rand::RandState, Integer};

use rhok;

const K: Range<u64> = 1..221;
const BITS: [u32; 4] = [64, 72, 84, 98];

fn random_prime(bits: u32, rng: &mut RandState) -> Integer {
    let mut p: Integer = Integer::random_bits(bits >> 1, rng).into();
    while p.is_probably_prime(48) == IsPrime::No {
        p = Integer::random_bits(bits >> 1, rng).into();
    }
    p
}

fn random_semiprime(bits: u32, rng: &mut RandState) -> Integer {
    random_prime(bits, rng) * random_prime(bits, rng)
}

fn bench_semiprime(c: &mut Criterion) {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(42u64 << 42));

    for bits in BITS {
        let mut group = c.benchmark_group(format!("semiprime-{}", bits));

        for k in K {
            group.bench_function(BenchmarkId::from_parameter(k), |b| {
                b.iter_batched(
                    || random_semiprime(bits, &mut rng),
                    |n| rhok::multi::pollard_rho(black_box(&n), black_box(k)),
                    criterion::BatchSize::SmallInput,
                )
            });
        }
    }
}

fn bench_general(c: &mut Criterion) {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(42u64 << 42));

    for bits in BITS {
        let mut group = c.benchmark_group(format!("general-{}", bits));

        for k in K {
            group.bench_function(BenchmarkId::from_parameter(k), |b| {
                b.iter_batched(
                    || {
                        random_prime(bits / 3, &mut rng)
                            * random_prime((2 * bits) / 3, &mut rng)
                    },
                    |n| rhok::multi::pollard_rho(black_box(&n), black_box(k)),
                    criterion::BatchSize::SmallInput,
                )
            });
        }
    }
}

fn bench_pow(c: &mut Criterion) {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(42u64 << 42));

    for bits in [4096] {
        let mut group = c.benchmark_group(format!("pow-{}", bits));
        group.warm_up_time(Duration::from_millis(100));

        for k in K {
            let _k = Integer::from(k << 1);
            group
                .bench_function(BenchmarkId::from_parameter(k), |b| {
                    b.iter_batched(
                        || {
                            (
                                Integer::from(Integer::random_bits(
                                    bits, &mut rng,
                                )),
                                Integer::from(Integer::random_bits(
                                    bits, &mut rng,
                                )),
                            )
                        },
                        |(a, b)| {
                            black_box(a).pow_mod(black_box(&_k), black_box(&b))
                        },
                        criterion::BatchSize::SmallInput,
                    )
                })
                .measurement_time(Duration::from_secs(1));
        }
    }
}

criterion_group!(
    name = factorize_benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500)).sample_size(50);
    targets = bench_general, bench_semiprime
);

criterion_group!(pow_benches, bench_pow);
criterion_main!(factorize_benches, pow_benches);
