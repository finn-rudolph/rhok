mod formula;
mod measurements;
mod pollard_rho;

use std::env;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Source {
    Formula,
    Measurements,
}

// Prints the running times / formula values for all possible assignments of
// k-values to a given number of machines. The k-values are maintained in `k`,
// `source` decides whether the formula or measurements are used.
fn iterate_k_values(
    k_min: u64,
    k_max: u64,
    source: Source,
    k: &mut Vec<u64>,
    j: usize,
    val: &mut Vec<f64>,
) {
    if j == k.len() {
        val.push(match source {
            Source::Measurements => measurements::measure_min_time(k),
            Source::Formula => formula::expected_time(k),
        });
        eprintln!("{:?}", k);
        return;
    }

    k[j] = if j > 0 { k[j - 1] } else { k_min };
    while k[j] <= k_max {
        iterate_k_values(k_min, k_max, source, k, j + 1, val);
        k[j] += 1;
    }
}

// TODO: maybe use higher order fn to not repeat the logic for iteration over k

fn print_values(
    k_min: u64,
    k_max: u64,
    k: &mut Vec<u64>,
    j: usize,
    val: &Vec<f64>,
    i: &mut usize,
    raw: bool,
) {
    if j == k.len() {
        if !raw {
            for k_j in k.iter() {
                print!("{:<5}", k_j);
            }
        }

        println!("{}", val[*i]);
        *i += 1;
        return;
    }

    k[j] = if j > 0 { k[j - 1] } else { k_min };
    while k[j] <= k_max {
        print_values(k_min, k_max, k, j + 1, val, i, raw);
        k[j] += 1;
    }
}

// CLI usage:
// cargo run --release -- [--formula | --measurements] [#machines]
//     [minimal k] [maximal k]
// If `--raw` is written after this, only the values without additional
// information are printed.

fn main() {
    let args: Vec<String> = env::args().collect();
    assert!(5 <= args.len() && args.len() <= 6);

    let source = match args[1].as_str() {
        "--formula" => Source::Formula,
        "--measurements" => Source::Measurements,
        _ => unreachable!(),
    };

    let machines: usize = args[2].parse().unwrap();
    let k_min: u64 = args[3].parse().unwrap();
    let k_max: u64 = args[4].parse().unwrap();

    assert!(machines > 0);
    assert!(k_min > 0);
    assert!(k_min <= k_max);

    let mut k = vec![0u64; machines];
    let mut values = Vec::new();

    iterate_k_values(k_min, k_max, source, &mut k, 0, &mut values);

    let raw = if args.len() == 6 {
        assert!(args[5].as_str() == "--raw");
        true
    } else {
        for i in 1..=machines {
            print!("{:<5}", format!("k{}", i));
        }
        println!("value");
        false
    };

    // Normalize values
    {
        let all_k_one = values[0];
        for v in values.iter_mut() {
            *v /= all_k_one;
        }
    }

    let mut i: usize = 0;
    k = vec![0u64; machines];
    print_values(k_min, k_max, &mut k, 0, &values, &mut i, raw);
}
