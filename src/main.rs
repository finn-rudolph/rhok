use std::{env, time::Duration};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rug::{Complete, Integer};

fn main() {
    rhok::stats::floyd_iteration_min_expectation_m2_gcd2();
}
