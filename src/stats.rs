use std::collections::VecDeque;

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

pub fn iota_stats() {
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

pub fn mu_lambda_nu_stats() {
    const A: usize = 1 << 16;
    const B: usize = (1 << 16) + (1 << 16);
    const K1: usize = 1;
    const K2: usize = 3;

    println!(
        "{:<20}{:<16}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
        "",
        "k",
        "gcd(2k, p - 1)",
        "nu mean",
        "nu std dev",
        "lambda mean",
        "lambda std dev",
        "mu mean",
        "mu std dev"
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
                "{:<20}{:<16}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
                "",
                k,
                gcd(2 * k as u64, p as u64 - 1),
                mean(&nu),
                standard_deviation(&nu),
                mean(&lambda),
                standard_deviation(&lambda),
                mean(&mu),
                standard_deviation(&mu)
            );
        }
    });
}

// Takes two sorted sequences of samples of random variables X, Y and returns
// the expected value of min(X, Y) assuming each (x, y)-pair is equally likely.
fn min_expectation_var2(x: &Vec<usize>, y: &Vec<usize>) -> f64 {
    let mut sum: usize = 0;
    let mut num_samples: usize = 0;

    for x_i in x {
        let (mut a, mut b) = (0, y.len());
        while a < b {
            let mid = (a + b) / 2;
            if y[mid] < *x_i {
                a = mid + 1;
            } else {
                b = mid;
            }
        }
        sum += (y.len() - a) * x_i;
        num_samples += a;
    }

    for y_i in y {
        let (mut a, mut b) = (0, x.len());
        while a < b {
            let mid = (a + b) / 2;
            if x[mid] <= *y_i {
                a = mid + 1;
            } else {
                b = mid;
            }
        }
        sum += (x.len() - a) * y_i;
        num_samples += a;
    }

    assert_eq!(num_samples, x.len() * y.len());
    (sum as f64) / ((x.len() * y.len()) as f64)
}

// The expected value of the number of steps until a collision occurs when
// running two machines in parallel and ignoring the higher costs of
// exponentiation when k != 1, but also requiring that the gcd(p - 1, 2k) = 2.
// Calculated over all possible assignments of k.
pub fn nu_min_expectation_m2_gcd2() {
    const A: usize = 1 << 12;
    const B: usize = (1 << 12) + (1 << 12);
    const K1: usize = 1;
    const K2: usize = 3;

    println!("{:<4} {:<4} {}", "k_1", "k_2", "E");

    (A..=B).for_each(|p| {
        if !miller_rabin(p as u64) {
            return;
        }

        println!("{}", p);
        let mtg = &Montgomery::new(p as u64);

        let nu: Vec<Vec<usize>> = (K1..=K2)
            .map(|k| {
                let (mu, lambda) = get_mu_lambda(p, k, &mtg);
                let mut nu: Vec<usize> =
                    mu.iter().zip(lambda.iter()).map(|(x, y)| x + y).collect();
                nu.sort();
                nu
            })
            .collect();

        for k1 in K1..=K2 {
            for k2 in k1..=K2 {
                if gcd(k1 as u64, (p as u64 - 1) / 2) == 1
                    && gcd(k2 as u64, (p as u64 - 1) / 2) == 1
                {
                    println!(
                        "{:<4} {:<4} {}",
                        k1,
                        k2,
                        min_expectation_var2(&nu[k1 - K1], &nu[k2 - K1])
                    );
                }
            }
        }
    });
}

// Same as above, but using the number of steps until a collision is detected
// by Floyd's algorithm instead.
pub fn floyd_iteration_min_expectation_m2_gcd2() {
    const A: usize = 2;
    const B: usize = 5;
    const K1: usize = 1;
    const K2: usize = 1;

    let floyd_min_m2: Vec<[[f64; K2 - K1 + 1]; K2 - K1 + 1]> = (A..=B)
        .into_par_iter()
        .map(|p| {
            let mut expectation = [[0.0; K2 - K1 + 1]; K2 - K1 + 1];

            if !miller_rabin(p as u64) {
                return expectation;
            }

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

                    assert!(floyd_iterations.iter().enumerate().all(
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
