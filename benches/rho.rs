use std::{ops::Range, time::Duration};

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rug::{integer::IsPrime, rand::RandState, Integer};

use rhok;

const K: Range<u32> = 2..1000;
const BITS: [u32; 4] = [64, 72, 84, 98];

fn random_prime(bits: u32, rng: &mut RandState) -> Integer {
    let mut p: Integer = Integer::random_bits(bits, rng).into();
    while p.is_probably_prime(48) == IsPrime::No {
        p = Integer::random_bits(bits, rng).into();
    }
    p
}

fn bench_semiprime(c: &mut Criterion) {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(42u64 << 42));

    for bits in BITS {
        let mut group = c.benchmark_group(format!("semiprime {}", bits));
        for k in K {
            group.bench_function(BenchmarkId::from_parameter(k), |b| {
                b.iter_batched(
                    || {
                        random_prime(bits >> 1, &mut rng)
                            * random_prime(bits >> 1, &mut rng)
                    },
                    |n| rhok::pollard_rho(black_box(&n), black_box(k)),
                    criterion::BatchSize::SmallInput,
                )
            });
        }
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_nanos(1))
        .sample_size(12);
    targets = bench_semiprime
);
criterion_main!(benches);
