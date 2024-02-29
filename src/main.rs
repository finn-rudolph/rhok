mod formula;
mod multi;
mod single;

use std::env;

use rand::RngCore;
use rug::Integer;

use crate::{multi::gen_test_numbers, single::miller_rabin};

const SAMPLES: usize = 1 << 10;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Source {
    Formula,
    Single,
    Multi,
}

fn random_prime(bits: u32, rng: &mut dyn RngCore) -> u64 {
    loop {
        let p = rng.next_u64() >> (64 - bits);
        if p != 2 && miller_rabin(p) {
            return p;
        }
    }
}

// Prints the running times / formula values for all possible assignments of
// k-values to a given number of machines. The k-values are maintained in `k`,
// `source` decides whether the formula or real measurements are used.
fn iterate_k_cartesian_product(
    k_min: u64,
    k_max: u64,
    source: Source,
    raw: bool,
    k: &mut Vec<u64>,
    j: usize,
    test_numbers: &Vec<Integer>,
) {
    if j == k.len() {
        if !raw {
            for k_j in k.iter() {
                print!("{:<5}", k_j);
            }
        }

        match source {
            Source::Single => println!("{}", single::measure(k)),
            Source::Multi => {
                let (val, outliers) = multi::measure(k, test_numbers);
                println!(
                    "{:<20}{:<20}{}",
                    (outliers as f64 / (SAMPLES + outliers) as f64) * 100.0,
                    outliers,
                    val
                )
            }
            Source::Formula => {
                println!("{}", formula::expected_time(k_min, k_max, k));
            }
        };

        return;
    }

    k[j] = if j > 0 { k[j - 1] } else { k_min };
    while k[j] <= k_max {
        iterate_k_cartesian_product(
            k_min,
            k_max,
            source,
            raw,
            k,
            j + 1,
            test_numbers,
        );
        k[j] += 1;
    }
}

// CLI usage: [--formula | --single | --multi] [#machines] [minimal k]
// [maximal k]
// If `--raw` is written after this, only the values without additional
// information are printed.

fn main() {
    // Only use half of the available threads for better measurement accuracy.
    rayon::ThreadPoolBuilder::new()
        .num_threads(10)
        .build_global()
        .unwrap();

    let args: Vec<String> = env::args().collect();
    assert!(5 <= args.len() && args.len() <= 6);

    let source = match args[1].as_str() {
        "--formula" => Source::Formula,
        "--single" => Source::Single,
        "--multi" => Source::Multi,
        _ => unreachable!(),
    };

    let machines: usize = args[2].parse().unwrap();
    let k_min: u64 = args[3].parse().unwrap();
    let k_max: u64 = args[4].parse().unwrap();

    if source != Source::Formula {
        println!("Number of (valid) samples per k-config: {}", SAMPLES);
    }

    let raw = if args.len() == 6 {
        assert!(args[5].as_str() == "--raw");
        true
    } else {
        for i in 1..=machines {
            print!("{:<5}", format!("k{}", i));
        }

        if source != Source::Formula {
            print!("{:<20}{:<20}", "% outliers", "# outliers");
        }

        println!("value");

        false
    };

    assert!(machines > 0);
    assert!(k_min > 0);
    assert!(k_min <= k_max);

    let mut k = vec![0u64; machines];

    iterate_k_cartesian_product(
        k_min,
        k_max,
        source,
        raw,
        &mut k,
        0,
        &gen_test_numbers(SAMPLES),
    );
}
