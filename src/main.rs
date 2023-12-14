use rand_xoshiro::rand_core::{RngCore, SeedableRng};

fn main() {
    let a = 5000;
    let b = 10000;

    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
    unsafe {
        for i in a..b {
            for k in 0..3 {
                let mut unused: Vec<usize> = (0..i).collect();
                for j in 0..i {
                    let next_index = rng.next_u64() as usize % unused.len();
                    rhok::stats::GRAPH[i - a][k][j] = unused[next_index];
                    unused.swap_remove(next_index);
                }
            }
        }
    }
    rhok::stats::nu_min_expectation_m2_gcd2::<1, 3>(a, b - 1);
}
