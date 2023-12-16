use std::collections::VecDeque;

use num_traits::AsPrimitive;
use rayon::prelude::*;

use crate::{
    fenwick_tree::FenwickTree,
    miller_rabin::miller_rabin,
    montgomery::Montgomery,
    single::{self, gcd},
};

fn f(x: usize, k: usize, mtg: &Montgomery) -> usize {
    let r = mtg.out_of_montgomery_space(mtg.add(
        mtg.pow(mtg.to_montgomery_space(x as u64), (k as u64) << 1),
        mtg.one(),
    ));
    assert!(r < 2 * mtg.n());
    (r - if r >= mtg.n() { mtg.n() } else { 0 }) as usize
}

fn create_histogram(x: &Vec<usize>) -> Vec<usize> {
    let mut histogram = vec![0; x.iter().max().unwrap() + 1];
    for v in x {
        histogram[*v] += 1;
    }
    histogram
}

pub fn histogram_to_distr(h: &[usize]) -> Vec<usize> {
    let mut x: Vec<usize> = Vec::new();
    for (i, f) in h.iter().enumerate() {
        for _ in 0..*f {
            x.push(i);
        }
    }
    x
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

    let mut variance = 0.0;
    for v in x {
        variance += (v.as_() - mean) * (v.as_() as f64 - mean);
    }

    (variance / x.len() as f64).sqrt()
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

fn get_cycle_nodes(
    start: usize,
    k: usize,
    mtg: &Montgomery,
    visited: &mut Vec<bool>,
) -> Vec<usize> {
    let mut cycle_nodes: Vec<usize> = Vec::new();

    let (mut x, mut y) = (start, f(start, k, &mtg));
    while x != y {
        x = f(x, k, &mtg);
        y = f(f(y, k, &mtg), k, &mtg);
    }

    while cycle_nodes.is_empty() || x != cycle_nodes[0] {
        cycle_nodes.push(x);
        visited[x] = true;
        x = f(x, k, &mtg);
    }

    cycle_nodes
}

fn get_connected_component(
    start: usize,
    k: usize,
    mtg: &Montgomery,
    pre: &Vec<Vec<usize>>,
    visited: &mut Vec<bool>,
) -> Vec<usize> {
    let mut q = VecDeque::new();
    let mut cc: Vec<usize> = Vec::new();
    q.push_back(start);
    visited[start] = true;

    while !q.is_empty() {
        let u = *q.front().unwrap();
        cc.push(u);
        q.pop_front();

        for v in &pre[u] {
            if !visited[*v] {
                visited[*v] = true;
                q.push_back(*v);
            }
        }

        let next = f(u, k, &mtg);
        if !visited[next] {
            visited[next] = true;
            q.push_back(next);
        }
    }

    cc
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
            let cycle_nodes = get_cycle_nodes(i, k, mtg, &mut visited);

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

// Returns (p, mu, lambda, nu), where in p are the primes from which the data
// was gathered.
fn get_mu_lambda_nu_mean<const K1: usize, const K2: usize>(
    a: usize,
    b: usize,
    g: usize,
    normalize: bool,
) -> [(Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1]
where
    [(); K2 - K1 + 1]:,
    [(Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1]: Default,
{
    let mln: Vec<(usize, [(f64, f64, f64); K2 - K1 + 1])> = (a..=b)
        .into_par_iter()
        .map(|p| {
            let mut mln = [(0.0, 0.0, 0.0); K2 - K1 + 1];

            if !miller_rabin(p as u64) {
                return (p, mln);
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

            (p, mln)
        })
        .collect();

    let mut result: [(Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1] =
        Default::default();

    for k in K1..=K2 {
        result[k - K1] = (
            mln.iter()
                .filter(|(_, x)| x[k - K1].0 != 0.0)
                .map(|(p, _)| *p)
                .collect::<Vec<usize>>(),
            mln.iter()
                .map(|(_, x)| x[k - K1].0)
                .filter(|x| *x != 0.0)
                .collect::<Vec<f64>>(),
            mln.iter()
                .map(|(_, x)| x[k - K1].1)
                .filter(|x| *x != 0.0)
                .collect::<Vec<f64>>(),
            mln.iter()
                .map(|(_, x)| x[k - K1].2)
                .filter(|x| *x != 0.0)
                .collect::<Vec<f64>>(),
        )
    }

    result
}

pub fn iota_individual(a: usize, b: usize, k_1: usize, k_2: usize) {
    println!(
        "{:<20}{:<16}{:<16}{:<20}{}",
        "", "k", "gcd(2k, p - 1)", "iota std dev", "iota histogram"
    );
    for p in a..=b {
        println!("p = {}", p,);
        let mtg = Montgomery::new(p as u64);

        for k in k_1..=k_2 {
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

pub fn mu_lambda_nu_individual(a: usize, b: usize, k_1: usize, k_2: usize) {
    // println!(
    //     "{:<12}{:<4}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{}\n",
    //     "",
    //     "k",
    //     "gcd(2k, p - 1)",
    //     "mu mean",
    //     "mu std dev",
    //     "lambda mean",
    //     "lambda std dev",
    //     "nu mean",
    //     "nu std dev",
    //     "nu histogram"
    // );

    let maxima: Vec<(usize, f64)> = (a..=b)
        .into_par_iter()
        .map(|p| {
            if !miller_rabin(p as u64) {
                return (p, 0.0);
            }

            // println!("{}", p);

            let mtg = Montgomery::new(p as u64);

            for k in k_1..=k_2 {
                let (mu, _) = get_mu_lambda(p, k, &mtg);
                let curr_max =
                    *mu.iter().max().unwrap() as f64 / (p as f64).sqrt();
                return (p, curr_max);

                // let nu: Vec<usize> =
                //     mu.iter().zip(lambda.iter()).map(|(x, y)| x + y).collect();

                // println!(
                //     "{:<12}{:<4}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20} {:?}",
                //     "",
                //     k,
                //     gcd(2 * k as u64, p as u64 - 1),
                //     mean(&mu),
                //     standard_deviation(&mu),
                //     mean(&lambda),
                //     standard_deviation(&lambda),
                //     mean(&nu),
                //     standard_deviation(&nu),
                //     create_histogram(&nu),
                // );
            }

            return (p, 0.0);
        })
        .filter(|x| x.1 != 0.0)
        .collect();

    let mut max_lambda_by_sqrt: f64 = 0.0;
    for (p, l) in maxima {
        if l > max_lambda_by_sqrt {
            println!("new max at p = {}, max = {}", p, l);
            max_lambda_by_sqrt = l;
        }
    }
}

// Calculates the mean and standard deviation of mu / lambda / nu for a fixed
// gcd over and for multiple k over primes in an interval.
pub fn mu_lambda_nu_summary<const K1: usize, const K2: usize>(
    a: usize,
    b: usize,
    k_gcd: usize,
    normalize: bool,
) where
    [(); K2 - K1 + 1]:,
    [(Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1]: Default,
{
    println!(
        "A = {}, B = {}, K1 = {}, K2 = {}\n\
        gcd(p - 1, 2k) = {}\nnormalization by sqrt p: {}\n",
        a, b, K1, K2, k_gcd, normalize
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

    let mln = get_mu_lambda_nu_mean::<K1, K2>(a, b, k_gcd, normalize);

    for k in K1..=K2 {
        let (_, mu, lambda, nu) = &mln[k - K1];
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

fn sorted_seq_intersection(x: &[usize], y: &[usize]) -> Vec<usize> {
    let mut z: Vec<usize> = Vec::new();
    let (mut i, mut j) = (x.iter().peekable(), y.iter().peekable());
    while i.peek().is_some() && j.peek().is_some() {
        if i.peek() < j.peek() {
            i.next().unwrap();
        } else if i.peek() > j.peek() {
            j.next().unwrap();
        } else {
            z.push(*i.next().unwrap());
            j.next().unwrap();
        }
    }
    z
}

// x is indexed by index and subseq_index is a subset of index. x[i] appears in
// the result if and only if index[i] is in subseq_index.
fn subsequence<T>(x: &[T], index: &[usize], subseq_index: &[usize]) -> Vec<T>
where
    T: Copy,
{
    let mut z: Vec<T> = Vec::new();
    let mut i = subseq_index.iter().peekable();
    for (v, j) in x.iter().zip(index.iter()) {
        if j == *i.peek().unwrap() {
            z.push(*v);
            i.next().unwrap();
            if i.peek().is_none() {
                break;
            }
        }
    }
    z
}

pub fn mu_lambda_nu_correlation<const K1: usize, const K2: usize>(
    a: usize,
    b: usize,
    k_gcd: usize,
    normalize: bool,
) where
    [(); K2 - K1 + 1]:,
    [(Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>); K2 - K1 + 1]: Default,
{
    println!(
        "A = {}, B = {}, K1 = {}, K2 = {}\n\
        gcd(p - 1, 2k) = {}\nnormalization by sqrt p: {}\n",
        a, b, K1, K2, k_gcd, normalize
    );

    println!(
        "{:<6}{:<6}{:<10}{:<24}{:<24}{:<24}\n",
        "k_1", "k_2", "#samples", "mu corr", "lambda corr", "nu corr",
    );

    let mln = get_mu_lambda_nu_mean::<K1, K2>(a, b, k_gcd, normalize);

    for k in K1..=K2 {
        for l in k + 1..=K2 {
            let (p1, mu1, lambda1, nu1) = &mln[k - K1];
            let (p2, mu2, lambda2, nu2) = &mln[l - K1];
            let q = sorted_seq_intersection(p1, p2);

            println!(
                "{:<6}{:<6}{:<10}{:<24}{:<24}{:<24}",
                k,
                l,
                q.len(),
                correlation(
                    &subsequence(mu1, p1, &q),
                    &subsequence(mu2, p2, &q)
                ),
                correlation(
                    &subsequence(lambda1, p1, &q),
                    &subsequence(lambda2, p2, &q)
                ),
                correlation(
                    &subsequence(nu1, p1, &q),
                    &subsequence(nu2, p2, &q)
                )
            )
        }
    }
}

pub fn cc_summary<const K1: usize, const K2: usize>(
    a: usize,
    b: usize,
    k_gcd: usize,
) where
    [(); K2 - K1 + 1]:,
    [(usize, Vec<usize>); K2 - K1 + 1]: Default,
{
    // println!(
    //     "A = {}, B = {}, K1 = {}, K2 = {}\n\
    //     gcd(p - 1, 2k) = {}\n",
    //     a, b, K1, K2, k_gcd
    // );

    // println!(
    //     "{:<6}{:<10}{:<24}{:<24}{:<24}{:<24}\n",
    //     "k",
    //     "#samples",
    //     "#cc mean",
    //     "#cc std dev",
    //     "cc size mean",
    //     "cc size std dev"
    // );

    let cc_num_size_distr: Vec<[(usize, Vec<usize>); K2 - K1 + 1]> = (a..=b)
        .collect::<Vec<usize>>()
        .into_par_iter()
        .map(|p| {
            let mut cc_num_sizes: [(usize, Vec<usize>); K2 - K1 + 1] =
                Default::default();
            if !miller_rabin(p as u64) {
                return cc_num_sizes;
            }

            let mtg = Montgomery::new(p as u64);

            for k in K1..=K2 {
                if gcd(p as u64 - 1, 2 * k as u64) as usize != k_gcd {
                    continue;
                }

                let pre = get_predecessors(p, k, &mtg);
                let mut visited = vec![false; p];
                let mut num_ccs: usize = 0;
                let mut cc_sizes: Vec<usize> = Vec::new();

                for x in 0..p {
                    if !visited[x] {
                        num_ccs += 1;
                        cc_sizes.push(
                            get_connected_component(
                                x,
                                k,
                                &mtg,
                                &pre,
                                &mut visited,
                            )
                            .len(),
                        );
                    }
                }

                cc_num_sizes[k - K1] = (num_ccs, cc_sizes);
            }

            cc_num_sizes
        })
        .filter(|x| x[0].0 != 0)
        .collect();

    for c in cc_num_size_distr {
        println!("{:?}", c);
    }

    // for k in K1..=K2 {
    //     let cc_num: Vec<usize> = cc_num_mean_size
    //         .iter()
    //         .map(|x| x[k - K1].0)
    //         .filter(|x| *x != 0)
    //         .collect();
    //     let cc_mean_size: Vec<f64> = cc_num_mean_size
    //         .iter()
    //         .map(|x| x[k - K1].1)
    //         .filter(|x| *x != 0.0)
    //         .collect();

    //     println!(
    //         "{:<6}{:<10}{:<24}{:<24}{:<24}{:<24}",
    //         k,
    //         cc_num.len(),
    //         mean(&cc_num),
    //         standard_deviation(&cc_num),
    //         mean(&cc_mean_size),
    //         standard_deviation(&cc_mean_size),
    //     );
    // }
}

// Takes two sorted sequences of samples of random variables X, Y and returns
// the expected value of min(X, Y) assuming each (x, y)-pair is equally likely.
fn min_expectation_var2<'a, T>(x: &[T], y: &[T]) -> f64
where
    T: std::cmp::PartialOrd + AsPrimitive<f64>,
{
    let mut sum: f64 = 0.0;

    let mut p: usize = 0; // pointer for the two pointer method
    for x_i in x {
        while p < y.len() && y[p] < *x_i {
            p += 1;
        }
        sum += <usize as AsPrimitive<f64>>::as_(y.len() - p) * x_i.as_()
    }

    p = 0;
    for y_i in y {
        while p < x.len() && x[p] <= *y_i {
            p += 1;
        }
        sum += <usize as AsPrimitive<f64>>::as_(x.len() - p) * y_i.as_();
    }

    (sum as f64) / ((x.len() * y.len()) as f64)
}

// The expected value of the number of steps until a collision occurs when
// running two machines in parallel and ignoring the higher costs of
// exponentiation when k != 1, but also requiring that the gcd(p - 1, 2k) = 2.
// Calculated over all possible assignments of k.
pub fn nu_min_expectation_m2_gcd2<const K1: usize, const K2: usize>(
    a: usize,
    b: usize,
) where
    [(); K2 - K1 + 1]:,
{
    // println!(
    //     "{:<12}{:<4} {:<4} {}",
    //     "", "k_1", "k_2", "E(min(nu_1, nu_2))"
    // );

    let nu_min_m2: Vec<[[f64; K2 - K1 + 1]; K2 - K1 + 1]> = (a..=b)
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
                            min_expectation_var2(&nu[k1 - K1], &nu[k2 - K1])
                                / (p as f64).sqrt();

                        // println!(
                        //     "{:<12}{:<4} {:<4} {}",
                        //     "",
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

pub fn disj(a: usize, b: usize, k: usize) {
    let proportion_disj: Vec<f64> = (a..=b)
        .into_par_iter()
        .map(|p| {
            if !miller_rabin(p as u64) {
                return -1.0;
            }

            let mtg = Montgomery::new(p as u64);
            let pre = &get_predecessors(p, k, &mtg);
            let mut num_non_disj: usize = 0;

            let (mu, lambda) = get_mu_lambda(p, k, &mtg);

            let mut tree: FenwickTree<isize> = FenwickTree::new(p);
            let mut visited = vec![false; p];

            for i in 0..p {
                if !visited[i] {
                    let cc =
                        get_connected_component(i, k, &mtg, pre, &mut visited);

                    // First find all mu-nu-pairs in the current connected
                    // component and sort them by nu. We'll sweep over nu.
                    let mut nu_mu: Vec<(usize, usize)> = cc
                        .iter()
                        .map(|node| (mu[*node] + lambda[*node], mu[*node]))
                        .collect();
                    nu_mu.sort();

                    for (nu, mu) in nu_mu.iter().rev() {
                        // Paths from different trees are non-disjoint when the
                        // other node hits the cycle before I collide.
                        num_non_disj += tree.prefix_sum(*nu) as usize;
                        tree.update(*mu, 1);
                    }

                    // We currently ignore that paths can also be non-disjoint
                    // when the other mu is > my nu but both are in the same
                    // tree.

                    // Clear the Fenwick Tree.
                    for (_, mu) in nu_mu {
                        tree.update(mu, -1);
                    }
                }
            }

            // Account for two possible orderings of each starting pair and
            // add p, since when both start at the same node, the paths are
            // also non-disjoint.
            ((2 * num_non_disj + p) as f64) / ((p * p) as f64)
        })
        .filter(|x| *x != -1.0)
        .collect();

    println!(
        "proportion of non-disjoint pairs of paths: {}",
        proportion_disj.iter().sum::<f64>() / proportion_disj.len() as f64
    );
}

// Same as above, but using the number of steps until a collision is detected
// by Floyd's algorithm instead.
pub fn floyd_iteration_min_expectation_m2_gcd2<
    const K1: usize,
    const K2: usize,
>(
    a: usize,
    b: usize,
) where
    [(); K2 - K1 + 1]:,
{
    // println!("{:<4} {:<4} {}", "k_1", "k_2", "E");

    let floyd_min_m2: Vec<[[f64; K2 - K1 + 1]; K2 - K1 + 1]> = (a..=b)
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

pub fn analyze_distributions<'a, T>(x: &'a [T], y: &'a [T])
where
    T: std::cmp::PartialOrd
        + std::cmp::Ord
        + AsPrimitive<f64>
        + std::iter::Sum<&'a T>,
{
    println!("mean x: {}", mean(x));
    println!("mean y: {}", mean(y));
    println!("std dev x: {}", standard_deviation(x));
    println!("std dev y: {}", standard_deviation(y));
    let mut _x: Vec<T> = x.iter().map(|a| a.clone()).collect();
    let mut _y: Vec<T> = y.iter().map(|a| a.clone()).collect();
    _x.sort();
    _y.sort();

    println!("min(x, y): {}", min_expectation_var2(&_x, &_y));
    println!("min(x, x): {}", min_expectation_var2(&_x, &_x));
    println!("min(y, y): {}", min_expectation_var2(&_y, &_y));
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_f() {
        for p in 3..10001 {
            if !miller_rabin(p) {
                continue;
            }

            let mtg = Montgomery::new(p);
            for i in 0..p {
                println!("{}, {}", p, i);
                assert_eq!(f(i as usize, 1, &mtg) as u64, (i * i + 1) % p);
            }
        }
    }
}
