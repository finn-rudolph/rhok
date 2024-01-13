mod miller_rabin;
mod montgomery;
mod pollard_rho;

use std::env;

#[derive(Clone, Copy)]
enum Source {
    Formula,
    Real,
}

fn iterate_k_cartesian_product(
    machines: usize,
    k_min: usize,
    k_max: usize,
    source: Source,
    k: &mut Vec<usize>,
    j: usize,
) {
    if j == machines {
        for k_j in k {
            print!("{:<5}", k_j);
        }

        println!(
            " | {}",
            match source {
                Source::Real => 0.0,
                Source::Formula => 0.0,
            }
        );
        return;
    }

    k[j] = if j > 0 { k[j - 1] } else { k_min };
    while k[j] <= k_max {
        iterate_k_cartesian_product(machines, k_min, k_max, source, k, j + 1);
        k[j] += 1;
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 5);

    let source = match args[1].as_str() {
        "--formula" => Source::Formula,
        "--rela" => Source::Real,
        _ => unreachable!(),
    };

    let machines: usize = args[2].parse().unwrap();
    let k_min: usize = args[3].parse().unwrap();
    let k_max: usize = args[4].parse().unwrap();

    let mut k = vec![0usize; machines];

    iterate_k_cartesian_product(machines, k_min, k_max, source, &mut k, 0);
}
