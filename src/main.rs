use std::{env, time::Duration};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rug::{Complete, Integer};

fn main() {
    rhok::stats::m2_nu_min_expectation();
}
