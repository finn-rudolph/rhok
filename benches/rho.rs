use std::{ops::Range, time::Duration};

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rug::{integer::IsPrime, rand::RandState, Integer};

use rhok;

const K: Range<u32> = 2..10;
const BITS: [u32; 4] = [64, 72, 84, 98];

fn random_semiprime(bits: u32, rng: &mut RandState) -> Integer {
    let mut p: Integer = Integer::random_bits(bits >> 1, rng).into();
    while p.is_probably_prime(48) == IsPrime::No {
        p = Integer::random_bits(bits >> 1, rng).into();
    }
    let mut q: Integer = Integer::random_bits(bits >> 1, rng).into();
    while q == p || q.is_probably_prime(48) == IsPrime::No {
        q = Integer::random_bits(bits >> 1, rng).into();
    }
    p * q
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

criterion_group!(
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(100))
        .sample_size(10);
    targets = bench_semiprime
);
criterion_main!(benches);
