use std::{collections::VecDeque, ops::AddAssign};

use num_traits::{AsPrimitive, Num};
use rayon::prelude::*;
use std::iter;

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

// Returns the cycle of the CC of `start` in order.
fn get_cycle_nodes(start: usize, k: usize, mtg: &Montgomery) -> Vec<usize> {
    let mut cycle_nodes: Vec<usize> = Vec::new();

    let (mut x, mut y) = (start, f(start, k, &mtg));
    while x != y {
        x = f(x, k, &mtg);
        y = f(f(y, k, &mtg), k, &mtg);
    }

    while cycle_nodes.is_empty() || x != cycle_nodes[0] {
        cycle_nodes.push(x);
        x = f(x, k, &mtg);
    }

    cycle_nodes
}

// Returns all nodes in the tree in the functional graph with root `cycle_node`
// in BFS order. (so `cycle_node` itself is at index 0, then all it's children,
// and so on)
fn get_tree(
    cycle_node: usize,
    cycle_pre: usize,
    pre: &Vec<Vec<usize>>,
) -> Vec<usize> {
    let mut tree: Vec<usize> = vec![cycle_node];
    tree.extend(&pre[cycle_node]);
    tree.remove(
        tree.iter()
            .skip(1)
            .position(|node| *node == cycle_pre)
            .unwrap()
            + 1,
    );

    let mut i = 1;
    while i < tree.len() {
        let u = tree[i];
        tree.extend(&pre[u]);
        i += 1;
    }

    tree
}

fn get_connected_component(
    start: usize,
    k: usize,
    mtg: &Montgomery,
    pre: &Vec<Vec<usize>>,
) -> Vec<usize> {
    let cycle_nodes = get_cycle_nodes(start, k, mtg);

    // Iterate over all cycle nodes and their predecessors and flatten the trees
    // in a single vector.
    cycle_nodes
        .iter()
        .zip(cycle_nodes[1..].iter().chain(iter::once(&cycle_nodes[0])))
        .flat_map(|(cycle_pre, cycle_node)| {
            get_tree(*cycle_node, *cycle_pre, pre)
        })
        .collect()
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
            let cycle_nodes = get_cycle_nodes(i, k, mtg);
            for node in &cycle_nodes {
                visited[*node] = true;
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
                        let cc = get_connected_component(x, k, &mtg, &pre);
                        cc_sizes.push(cc.len());
                        for node in cc {
                            visited[node] = true;
                        }
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

// Takes two sorted sequences of numbers x, y and returns the sum of all minima
// of pairs (a, b), where a is in x and b is in y.
fn min_sum_var2<'a, T>(x: &[T], y: &[T]) -> T
where
    T: std::cmp::PartialOrd + Num + AddAssign<T> + Copy + 'static,
    usize: AsPrimitive<T>,
{
    let mut sum: T = T::zero();

    let mut p: usize = 0; // pointer for the two pointer method
    for x_i in x {
        while p < y.len() && y[p] < *x_i {
            p += 1;
        }
        sum += AsPrimitive::<T>::as_(y.len() - p) * *x_i
    }

    p = 0;
    for y_i in y {
        while p < x.len() && x[p] <= *y_i {
            p += 1;
        }
        sum += (x.len() - p).as_() * *y_i;
    }

    sum
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
                            min_sum_var2(&nu[k1 - K1], &nu[k2 - K1]) as f64
                                / ((p * p) as f64 * (p as f64).sqrt());

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

pub fn count_disjoint_paths(a: usize, b: usize, k: usize) {
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
            let nu: Vec<usize> = mu
                .iter()
                .zip(lambda.iter())
                .map(|(mu, lambda)| *mu + *lambda)
                .collect();

            let mut ftree: FenwickTree<isize> = FenwickTree::new(p + 1);
            let mut visited = vec![false; p];
            let mut is_cycle_node = vec![false; p];
            let mut num_desc_with_dis: Vec<Vec<usize>> = vec![Vec::new(); p];

            for i in 0..p {
                if !visited[i] {
                    let cc = get_connected_component(i, k, &mtg, pre);
                    for node in &cc {
                        visited[*node] = true;
                    }

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
                        num_non_disj += ftree.prefix_sum(*nu) as usize;
                        ftree.update(*mu, 1);
                    }

                    // Clear the Fenwick Tree.
                    for (_, mu) in nu_mu {
                        ftree.update(mu, -1);
                    }

                    // Now count the number of nodes in the same tree whose mu
                    // is > our nu, but who still reach our path.
                    let cycle_nodes = get_cycle_nodes(i, k, &mtg);
                    for node in &cycle_nodes {
                        is_cycle_node[*node] = true;
                    }
                    for (cycle_pre, cycle_node) in cycle_nodes.iter().zip(
                        cycle_nodes[1..]
                            .iter()
                            .chain(iter::once(&cycle_nodes[0])),
                    ) {
                        let tree = get_tree(*cycle_node, *cycle_pre, pre);

                        // Sweep over nu descending and keep in each node, how
                        // many descendants with a particular distance it has.
                        for node in tree.iter().rev() {
                            let mut u = *node;

                            // Go over all ancestors, collect the number of
                            // nodes that will hit the path of `node` and
                            // update the counts.
                            let mut i = 0;
                            while !is_cycle_node[u] {
                                // Only count the nodes that _just_ manage to
                                // reach that ancestor. (other nodes will be
                                // counted by higher ancestors already)
                                if num_desc_with_dis[u].len() > nu[u] {
                                    num_non_disj += num_desc_with_dis[u][nu[u]];
                                }
                                if num_desc_with_dis[u].len() <= i {
                                    num_desc_with_dis[u].resize(i + 1, 0);
                                }
                                num_desc_with_dis[u][i] += 1;

                                u = f(u, k, &mtg);
                                i += 1;
                            }
                        }

                        for node in tree {
                            num_desc_with_dis[node] = Vec::new();
                        }
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

// Computes the mean of min(nu_1, nu_2) for disjoint and non-disjoint paths
// separately. Runs in expected O(p sqrt p) per prime.
pub fn min_nu_disjoint_paths(a: usize, b: usize, k: usize) {
    // mean of disjoint, mean of non-disjoint, proportion of non-disjoint
    let disj_non_disj: Vec<(f64, f64, f64)> = (a..=b)
        .into_par_iter()
        .map(|p| {
            if !miller_rabin(p as u64) {
                return (-1.0, 0.0, 0.0);
            }

            let mtg = Montgomery::new(p as u64);
            let pre = &get_predecessors(p, k, &mtg);

            let mut num_non_disj: usize = 0;
            let mut nu_sum_non_disj: usize = 0;

            let (mu, lambda) = get_mu_lambda(p, k, &mtg);
            let mut nu: Vec<usize> = mu
                .iter()
                .zip(lambda.iter())
                .map(|(mu, lambda)| *mu + *lambda)
                .collect();

            let mut ftree: FenwickTree<isize> = FenwickTree::new(p + 1);
            let mut visited = vec![false; p];
            let mut is_cycle_node = vec![false; p];
            let mut num_desc_with_dis: Vec<Vec<usize>> = vec![Vec::new(); p];

            for i in 0..p {
                if !visited[i] {
                    let cc = get_connected_component(i, k, &mtg, pre);
                    for node in &cc {
                        visited[*node] = true;
                    }

                    // First find all mu-nu-pairs in the current connected
                    // component and sort them by nu. We'll sweep over nu.
                    let mut nu_mu: Vec<(usize, usize)> = cc
                        .iter()
                        .map(|node| (mu[*node] + lambda[*node], mu[*node]))
                        .collect();
                    nu_mu.sort();
                    assert!(nu_mu
                        .iter()
                        .all(|(nu, mu)| *nu - *mu == nu_mu[0].0 - nu_mu[0].1));

                    for (nu, mu) in nu_mu.iter().rev() {
                        // Paths are non-disjoint when the other node hits the
                        // cycle before I collide.
                        let num_cycle_coll = ftree.prefix_sum(*nu) as usize;
                        num_non_disj += num_cycle_coll;
                        nu_sum_non_disj += num_cycle_coll * nu;
                        ftree.update(*mu, 1);
                    }

                    // Clear the Fenwick Tree.
                    for (_, mu) in nu_mu {
                        ftree.update(mu, -1);
                    }

                    // Now count the number of nodes in the same tree whose mu
                    // is > our nu, but who still reach our path.
                    let cycle_nodes = get_cycle_nodes(i, k, &mtg);
                    for node in &cycle_nodes {
                        is_cycle_node[*node] = true;
                    }
                    for (cycle_pre, cycle_node) in cycle_nodes.iter().zip(
                        cycle_nodes[1..]
                            .iter()
                            .chain(iter::once(&cycle_nodes[0])),
                    ) {
                        let tree = get_tree(*cycle_node, *cycle_pre, pre);

                        // Sweep over nu descending and keep in each node, how
                        // many descendants with a particular distance it has.
                        for node in tree.iter().rev() {
                            let mut u = *node;

                            // Go over all ancestors, collect the number of
                            // nodes that will hit the path of `node` and
                            // update the counts.
                            let mut i = 0;
                            while !is_cycle_node[u] {
                                // Only count the nodes that _just_ manage to
                                // reach that ancestor. (other nodes will be
                                // counted by higher ancestors already)
                                if num_desc_with_dis[u].len() > nu[*node] {
                                    num_non_disj +=
                                        num_desc_with_dis[u][nu[*node]];
                                    nu_sum_non_disj +=
                                        nu[u] * num_desc_with_dis[u][nu[*node]];
                                }
                                if num_desc_with_dis[u].len() <= i {
                                    num_desc_with_dis[u].resize(i + 1, 0);
                                }
                                num_desc_with_dis[u][i] += 1;

                                u = f(u, k, &mtg);
                                i += 1;
                            }
                        }

                        for node in tree {
                            num_desc_with_dis[node] = Vec::new();
                        }
                    }
                }
            }

            // Account for two possible orderings of each starting pair and
            // add p, since when both start at the same node, the paths are
            // also non-disjoint.
            num_non_disj = 2 * num_non_disj + p;
            assert!(num_non_disj <= p * p);
            nu_sum_non_disj = 2 * nu_sum_non_disj + nu.iter().sum::<usize>();

            nu.sort();
            let total_nu_sum = min_sum_var2(&nu, &nu);
            assert!(nu_sum_non_disj <= total_nu_sum);

            // If there are no disjoint paths, return -1.0 so the sample gets
            // filtered out.
            (
                if num_non_disj != p * p {
                    (((total_nu_sum - nu_sum_non_disj) as f64)
                        / (p * p - num_non_disj) as f64)
                        / (p as f64).sqrt()
                } else {
                    -1.0
                },
                nu_sum_non_disj as f64
                    / (num_non_disj as f64 * (p as f64).sqrt()),
                num_non_disj as f64 / ((p * p) as f64),
            )
        })
        .filter(|x| x.0 != -1.0)
        .collect();

    let num_samples = disj_non_disj.len() as f64;
    let (mean_disj, mean_non_disj, disj_ratio) = disj_non_disj
        .into_iter()
        .reduce(|acc, elem| (acc.0 + elem.0, acc.1 + elem.1, acc.2 + elem.2))
        .unwrap();

    println!(
        "mean disjoint = {}\nmean non-disjoint = {}\ndisjoint ratio = {}\n\
        #samples = {}",
        mean_disj / num_samples,
        mean_non_disj / num_samples,
        disj_ratio / num_samples,
        num_samples
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
                        expectation[k1 - K1][k2 - K1] = min_sum_var2(
                            &floyd_iterations[k1 - K1],
                            &floyd_iterations[k2 - K1],
                        )
                            as f64
                            / (p * p) as f64;

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
