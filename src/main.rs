use std::{env, time::Duration};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rug::{Complete, Integer};

fn main() {
    rhok::stats::mu_lambda_nu_summary();
}
