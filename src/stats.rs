use std::{collections::VecDeque, ops::Sub};

use num_traits::AsPrimitive;
use rayon::prelude::*;

use crate::{
    miller_rabin::miller_rabin,
    montgomery::Montgomery,
    single::{self, gcd},
};

fn f(x: usize, k: usize, mtg: &Montgomery) -> usize {
    mtg.strict(mtg.pow(x as u64, (k as u64) << 1) + 1) as usize
}

fn create_histogram(x: &Vec<usize>) -> Vec<usize> {
    let mut histogram = vec![0; x.iter().max().unwrap() + 1];
    for v in x {
        histogram[*v] += 1;
    }
    histogram
}

fn mean<'a, T>(x: &'a [T]) -> f64
where
    T: std::iter::Sum<&'a T> + AsPrimitive<f64>,
{
    T::from(x.iter().sum::<T>()).as_() / x.len() as f64
}

fn standard_deviation<'a, T>(x: &'a [T]) -> f64
where
    T: std::iter::Sum<&'a T> + AsPrimitive<f64>,
{
    let mean = mean(x);

    let mut sigma = 0.0;
    for v in x {
        sigma += (v.as_() - mean) * (v.as_() as f64 - mean);
    }
    (sigma / x.len() as f64).sqrt()
}

// The pearson correlation for two sequences x_i, y_i of samples of random
// variables, where each pair (x_i, y_i) is equally likely and (x_i, y_j) for
// i != j has probability 0.
fn correlation<'a, T>(x: &'a [T], y: &'a [T]) -> f64
where
    T: std::iter::Sum<&'a T> + AsPrimitive<f64>,
{
    assert_eq!(x.len(), y.len());
    let (mean_x, mean_y) = (mean(x), mean(y));
    let mut cov = 0.0;
    for (x_i, y_i) in x.iter().zip(y.iter()) {
        cov += (x_i.as_() - mean_x) * (y_i.as_() - mean_y);
    }
    cov / (x.len() as f64 * standard_deviation(x) * standard_deviation(y))
}

fn get_predecessors(p: usize, k: usize, mtg: &Montgomery) -> Vec<Vec<usize>> {
    let mut pre: Vec<Vec<usize>> = vec![Vec::new(); p];

    for i in 0..p {
        let next = f(i, k, mtg);
        pre[next].push(i);
    }

    pre
}

// iota(u) = # predecessors of u that have in-degree > 1 (i.e lie in the
// interesting subgraph)
fn get_iota(p: usize, k: usize, mtg: &Montgomery) -> Vec<usize> {
    let mut iota = vec![0; p];
    let pre = get_predecessors(p, k, &mtg);

    for i in 0..p {
        // there's some weird in-degree 1 node, we wanna ignore this
        if pre[i].len() > 1 {
            for u in &pre[i] {
                if pre[*u].len() > 1 {
                    iota[*u] += 1;
                }
            }
        }
    }

    iota
}

fn get_mu_lambda(
    p: usize,
    k: usize,
    mtg: &Montgomery,
) -> (Vec<usize>, Vec<usize>) {
    let (mut mu, mut lambda) = (vec![0usize; p], vec![0usize; p]);

    let pre = get_predecessors(p, k, &mtg);
    let mut visited: Vec<bool> = vec![false; p];

    for i in 0..p {
        if !visited[i] {
            let mut cycle_nodes: Vec<usize> = Vec::new();

            {
                // mark cycle nodes
                let (mut x, mut y) = (i, f(i, k, &mtg));
                while x != y {
                    x = f(x, k, &mtg);
                    y = f(f(y, k, &mtg), k, &mtg);
                }

                while cycle_nodes.is_empty() || x != cycle_nodes[0] {
                    cycle_nodes.push(x);
                    visited[x] = true;
                    x = f(x, k, &mtg);
                }
            }

            for x in &cycle_nodes {
                let mut q: VecDeque<(usize, usize)> = VecDeque::new();
                q.push_back((*x, 0));

                while !q.is_empty() {
                    let (u, dis) = *q.front().unwrap();
                    q.pop_front();
                    visited[u] = true;
                    mu[u] = dis;
                    lambda[u] = cycle_nodes.len();

                    for v in &pre[u] {
                        if !visited[*v] {
                            q.push_back((*v, dis + 1));
                        }
                    }
                }
            }
        }
    }

    (mu, lambda)
}

fn get_mu_lambda_nu_mean<const K1: usize, const K2: usize>(
    a: usize,
    b: usize,
    g: usize,
    normalize: bool,
) -> [(Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1]
where
    [(); K2 - K1 + 1]:,
    [(Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1]: Default,
{
    let mln: Vec<[(f64, f64, f64); K2 - K1 + 1]> = (a..=b)
        .into_par_iter()
        .map(|p| {
            let mut mln = [(0.0, 0.0, 0.0); K2 - K1 + 1];

            if !miller_rabin(p as u64) {
                return mln;
            }

            let mtg = Montgomery::new(p as u64);

            for k in K1..=K2 {
                if gcd(2 * k as u64, p as u64 - 1) as usize == g {
                    let (mu, lambda) = get_mu_lambda(p, k, &mtg);
                    let nu: Vec<usize> = mu
                        .iter()
                        .zip(lambda.iter())
                        .map(|(x, y)| x + y)
                        .collect();

                    let normalization_factor = if normalize {
                        1.0 / (p as f64).sqrt()
                    } else {
                        1.0
                    };

                    mln[k - K1] = (
                        mean(&mu) * normalization_factor,
                        mean(&lambda) * normalization_factor,
                        mean(&nu) * normalization_factor,
                    );
                }
            }

            mln
        })
        .collect();

    let mut result: [(Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1] =
        Default::default();

    for k in K1..=K2 {
        result[k - K1] = (
            mln.iter()
                .map(|x| x[k - K1].0)
                .filter(|x| *x != 0.0)
                .collect::<Vec<f64>>(),
            mln.iter()
                .map(|x| x[k - K1].1)
                .filter(|x| *x != 0.0)
                .collect::<Vec<f64>>(),
            mln.iter()
                .map(|x| x[k - K1].2)
                .filter(|x| *x != 0.0)
                .collect::<Vec<f64>>(),
        )
    }

    result
}

pub fn iota_individual() {
    const A: usize = 1 << 16;
    const B: usize = (1 << 16) + (1 << 12);
    const K1: usize = 1;
    const K2: usize = 12;

    println!(
        "{:<20}{:<16}{:<16}{:<20}{}",
        "", "k", "gcd(2k, p - 1)", "iota std dev", "iota histogram"
    );
    for p in A..=B {
        println!("p = {}", p,);
        let mtg = Montgomery::new(p as u64);

        for k in K1..=K2 {
            let iota = get_iota(p, k, &mtg);
            let iota_hist = create_histogram(&iota);

            let g = gcd(p as u64 - 1, 2 * k as u64);
            println!(
                "{:<20}{:<16}{:<16}{:<20} {:<?}",
                "",
                k,
                g,
                standard_deviation(&iota),
                &iota_hist,
            );

            if g == 2 {
                if p & 7 == 1 {
                    assert!(iota_hist[0] == (p - 1) >> 3);
                    assert!(iota_hist[1] == (p + 3) >> 2);
                    assert!(iota_hist[2] == (p - 9) >> 3);
                } else if p & 7 == 3 {
                    assert!(iota_hist[0] == (p - 3) >> 3);
                    assert!(iota_hist[1] == (p + 1) >> 2);
                    assert!(iota_hist[2] == (p - 3) >> 3);
                } else if p & 7 == 5 {
                    assert!(iota_hist[0] == (p + 3) >> 3);
                    assert!(iota_hist[1] == (p - 1) >> 2);
                    assert!(iota_hist[2] == (p - 5) >> 3);
                } else if p & 7 == 7 {
                    assert!(iota_hist[0] == (p + 1) >> 3);
                    assert!(iota_hist[1] == (p - 3) >> 2);
                    assert!(iota_hist[2] == (p + 1) >> 3);
                }
            }
        }
    }
}

pub fn mu_lambda_nu_individual() {
    const A: usize = 1 << 16;
    const B: usize = 1 << 17;
    const K1: usize = 1;
    const K2: usize = 3;

    println!(
        "{:<12}{:<4}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}\n",
        "",
        "k",
        "gcd(2k, p - 1)",
        "mu mean",
        "mu std dev",
        "lambda mean",
        "lambda std dev",
        "nu mean",
        "nu std dev"
    );

    (A..=B).for_each(|p| {
        if !miller_rabin(p as u64) {
            return;
        }

        println!("{}", p);

        let mtg = Montgomery::new(p as u64);

        for k in K1..=K2 {
            let (mu, lambda) = get_mu_lambda(p, k, &mtg);
            let nu: Vec<usize> =
                mu.iter().zip(lambda.iter()).map(|(x, y)| x + y).collect();

            println!(
                "{:<12}{:<4}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
                "",
                k,
                gcd(2 * k as u64, p as u64 - 1),
                mean(&mu),
                standard_deviation(&mu),
                mean(&lambda),
                standard_deviation(&lambda),
                mean(&nu),
                standard_deviation(&nu)
            );
        }
    });
}

// Calculates the mean and standard deviation of mu / lambda / nu for a fixed
// gcd over and for multiple k over primes in an interval.
pub fn mu_lambda_nu_summary() {
    const A: usize = 1 << 13;
    const B: usize = 1 << 14;
    const K1: usize = 1;
    const K2: usize = 30;
    const GCD: usize = 2;
    const NORMALIZE: bool = true;

    println!(
        "A = {}, B = {}, K1 = {}, K2 = {}\n\
        gcd(p - 1, 2k) = {}\nnormalization by sqrt p: {}\n",
        A, B, K1, K2, GCD, NORMALIZE
    );

    println!(
        "{:<4}{:<10}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}\n",
        "k",
        "#samples",
        "mu mean",
        "mu std dev",
        "lambda mean",
        "lambda std dev",
        "nu mean",
        "nu std dev",
    );

    let mln = get_mu_lambda_nu_mean::<K1, K2>(A, B, GCD, NORMALIZE);

    for k in K1..=K2 {
        let (mu, lambda, nu) = &mln[k - K1];
        assert_eq!(mu.len(), lambda.len());
        assert_eq!(lambda.len(), nu.len());

        println!(
            "{:<4}{:<10}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
            k,
            mu.len(),
            mean(&mu),
            standard_deviation(&mu),
            mean(&lambda),
            standard_deviation(&lambda),
            mean(&nu),
            standard_deviation(&nu)
        );
    }
}

// Takes two sorted sequences of samples of random variables X, Y and returns
// the expected value of min(X, Y) assuming each (x, y)-pair is equally likely.
fn min_expectation_var2(x: &Vec<usize>, y: &Vec<usize>) -> f64 {
    let mut sum: usize = 0;

    let mut p = 0; // pointer for the two pointer method
    for x_i in x {
        while p < y.len() && y[p] < *x_i {
            p += 1;
        }
        sum += (y.len() - p) * x_i;
    }

    p = 0;
    for y_i in y {
        while p < x.len() && x[p] <= *y_i {
            p += 1;
        }
        sum += (x.len() - p) * y_i;
    }

    (sum as f64) / ((x.len() * y.len()) as f64)
}

// The expected value of the number of steps until a collision occurs when
// running two machines in parallel and ignoring the higher costs of
// exponentiation when k != 1, but also requiring that the gcd(p - 1, 2k) = 2.
// Calculated over all possible assignments of k.
pub fn nu_min_expectation_m2_gcd2() {
    const A: usize = 1 << 15;
    const B: usize = 1 << 16;
    const K1: usize = 1;
    const K2: usize = 3;

    // println!("{:<4} {:<4} {}", "k_1", "k_2", "E");

    let nu_min_m2: Vec<[[f64; K2 - K1 + 1]; K2 - K1 + 1]> = (A..=B)
        .into_par_iter()
        .map(|p| {
            let mut expectation = [[0.0; K2 - K1 + 1]; K2 - K1 + 1];

            if !miller_rabin(p as u64) {
                return expectation;
            }

            // println!("{}", p);
            let mtg = &Montgomery::new(p as u64);

            let nu: Vec<Vec<usize>> = (K1..=K2)
                .map(|k| {
                    let (mu, lambda) = get_mu_lambda(p, k, &mtg);
                    let mut nu: Vec<usize> = mu
                        .iter()
                        .zip(lambda.iter())
                        .map(|(x, y)| x + y)
                        .collect();
                    nu.sort();
                    nu
                })
                .collect();

            for k1 in K1..=K2 {
                for k2 in k1..=K2 {
                    if gcd(k1 as u64, (p as u64 - 1) / 2) == 1
                        && gcd(k2 as u64, (p as u64 - 1) / 2) == 1
                    {
                        expectation[k1 - K1][k2 - K1] =
                            min_expectation_var2(&nu[k1 - K1], &nu[k2 - K1]);

                        // println!(
                        //     "{:<4} {:<4} {}",
                        //     k1,
                        //     k2,
                        //     expectation[k1 - K1][k2 - K1]
                        // );
                    }
                }
            }

            expectation
        })
        .filter(|x| !x.iter().all(|y| y.iter().all(|z| *z == 0.0)))
        .collect();

    println!(
        "{:<4} {:<4} {:<10} {:<20} {:<20}",
        "k_1", "k_2", "#samples", "mean nu", "std dev nu"
    );

    for k1 in K1..=K2 {
        for k2 in k1..=K2 {
            let u: Vec<f64> = nu_min_m2
                .iter()
                .map(|x| x[k1 - K1][k2 - K1])
                .filter(|x| *x != 0.0)
                .collect();

            println!(
                "{:<4} {:<4} {:<10} {:<20} {:<20}",
                k1,
                k2,
                u.len(),
                mean(&u),
                standard_deviation(&u)
            );
        }
    }
}

// Same as above, but using the number of steps until a collision is detected
// by Floyd's algorithm instead.
pub fn floyd_iteration_min_expectation_m2_gcd2() {
    const A: usize = 1 << 16;
    const B: usize = 1 << 17;
    const K1: usize = 1;
    const K2: usize = 3;

    // println!("{:<4} {:<4} {}", "k_1", "k_2", "E");

    let floyd_min_m2: Vec<[[f64; K2 - K1 + 1]; K2 - K1 + 1]> = (A..=B)
        .into_par_iter()
        .map(|p| {
            let mut expectation = [[0.0; K2 - K1 + 1]; K2 - K1 + 1];

            if !miller_rabin(p as u64) {
                return expectation;
            }

            // println!("{}", p);
            let mtg = Montgomery::new(p as u64);

            let floyd_iterations: Vec<Vec<usize>> = (K1..=K2)
                .map(|k| {
                    // A collision is detected in the c-th iteration, where c is
                    // minimal such that c = 2c mod lambda. But then it is easy to
                    // see that c must be ceil(mu / lambda) * lambda.

                    let (mu, lambda) = get_mu_lambda(p, k, &mtg);
                    let mut floyd_iterations = mu
                        .iter()
                        .zip(lambda.iter())
                        .map(|(mu, lambda)| {
                            (*lambda).max(((mu + lambda - 1) / lambda) * lambda)
                        })
                        .collect::<Vec<usize>>();

                    debug_assert!(floyd_iterations.iter().enumerate().all(
                        |(start, iterations)| {
                            *iterations
                                == single::pollard_slow_iteration_count(
                                    p as u64,
                                    k as u64,
                                    start as u64,
                                )
                        },
                    ));

                    floyd_iterations.sort();
                    floyd_iterations
                })
                .collect();

            for k1 in K1..=K2 {
                for k2 in k1..=K2 {
                    if gcd(k1 as u64, (p as u64 - 1) / 2) == 1
                        && gcd(k2 as u64, (p as u64 - 1) / 2) == 1
                    {
                        expectation[k1 - K1][k2 - K1] = min_expectation_var2(
                            &floyd_iterations[k1 - K1],
                            &floyd_iterations[k2 - K1],
                        );

                        // println!(
                        //     "{:<4} {:<4} {}",
                        //     k1,
                        //     k2,
                        //     expectation[k1 - K1][k2 - K1]
                        // );
                    }
                }
            }

            expectation
        })
        .filter(|x| !x.iter().all(|y| y.iter().all(|z| *z == 0.0)))
        .collect();

    println!(
        "{:<4} {:<4} {:<10} {:<20} {:<20}",
        "k_1", "k_2", "#samples", "mean floyd iter", "std dev floyd iter"
    );

    for k1 in K1..=K2 {
        for k2 in k1..=K2 {
            let u: Vec<f64> = floyd_min_m2
                .iter()
                .map(|x| x[k1 - K1][k2 - K1])
                .filter(|x| *x != 0.0)
                .collect();

            println!(
                "{:<4} {:<4} {:<10} {:<20} {:<20}",
                k1,
                k2,
                u.len(),
                mean(&u),
                standard_deviation(&u)
            );
        }
    }
}
