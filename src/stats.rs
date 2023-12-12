use std::collections::VecDeque;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{miller_rabin::miller_rabin, montgomery::Montgomery, single::gcd};

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

fn mean(x: &Vec<usize>) -> f64 {
    x.iter().sum::<usize>() as f64 / x.len() as f64
}

fn standard_deviation(x: &Vec<usize>) -> f64 {
    let mean = mean(x);
    let mut sigma = 0.0;
    for v in x {
        sigma += (*v as f64 - mean) * (*v as f64 - mean);
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

pub fn iota() {
    const A: usize = 1 << 16;
    const B: usize = (1 << 16) + (1 << 12);
    const K1: usize = 1;
    const K2: usize = 12;

    println!(
        "{:<20}{:<16}{:<16}{:<20}{}",
        "", "k", "gcd(2k, p - 1)", "std dev(iota)", "avg of squares (iota)"
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

// calculates average collision length (nu) for various k
pub fn mu_lambda() {
    const A: usize = 1 << 16;
    const B: usize = (1 << 16) + (1 << 12);
    const K1: usize = 1;
    const K2: usize = 12;

    println!(
        "{:<20}{:<16}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
        "",
        "k",
        "gcd(k, p - 1)",
        "nu mean",
        "nu std dev",
        "cycle mean",
        "cycle std dev",
        "tail mean",
        "tail std dev"
    );

    (A..=B).into_iter().for_each(|p| {
        let mut nu = [0.0; K2 + 1];
        let mut gcds = [0usize; K2 + 1];

        if !miller_rabin(p as u64) {
            return;
        }

        let mtg = Montgomery::new(p as u64);
        let mut collision_len: Vec<Vec<usize>> = vec![vec![0; p]; K2 + 1];
        let mut cycle_len: Vec<usize> = vec![0; p];
        let mut histogram: Vec<Vec<usize>> = vec![vec![0; p]; K2 + 1];
        let mut cycle_avg = [0f64; K2 + 1];

        println!("{}", p);

        for k in K1..=K2 {
            gcds[k] = gcd(2 * k as u64, p as u64 - 1) as usize;

            let mut total_collision_len: usize = 0;
            let mut total_cycle_len: usize = 0;

            let mean_nu = (total_collision_len as f64) / (p as f64);
            nu[k] = mean_nu;
            let mean_cycle_len = (total_cycle_len as f64) / (p as f64);
            let mean_tail_len =
                ((total_collision_len - total_cycle_len) as f64) / (p as f64);

            let mut std_dev_nu = 0.0;
            let mut std_dev_cycle_len = 0.0;
            let mut std_dev_tail_len = 0.0;
            for i in 0..p {
                std_dev_nu += (collision_len[k][i] as f64 - mean_nu)
                    * (collision_len[k][i] as f64 - mean_nu);
                std_dev_cycle_len += (cycle_len[i] as f64 - mean_cycle_len)
                    * (cycle_len[i] as f64 - mean_cycle_len);
                std_dev_tail_len += ((collision_len[k][i] - cycle_len[i])
                    as f64
                    - mean_tail_len)
                    * ((collision_len[k][i] - cycle_len[i]) as f64
                        - mean_tail_len);
            }

            let g = gcd(2 * k as u64, p as u64 - 1);
            if g == 2 {
                cycle_avg[k] = mean_cycle_len;
            }

            std_dev_nu = (std_dev_nu / p as f64).sqrt();
            std_dev_cycle_len = (std_dev_cycle_len / p as f64).sqrt();
            std_dev_tail_len = (std_dev_tail_len / p as f64).sqrt();

            println!(
                "{:<20}{:<16}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
                "",
                k,
                g,
                mean_nu,
                std_dev_nu,
                mean_cycle_len,
                std_dev_cycle_len,
                mean_tail_len,
                std_dev_tail_len
            );

            collision_len[k].sort();
        }

        for k1 in K1..=K2 {
            for k2 in k1..=K2 {
                if gcd(k1 as u64, (p as u64 - 1) / 2) == 1
                    && gcd(k2 as u64, (p as u64 - 1) / 2) == 1
                {
                    let mut sum: usize = 0;
                    let mut num_samples: usize = 0;

                    for x in &collision_len[k1] {
                        let (mut a, mut b) = (0, p);
                        while a < b {
                            let mid = (a + b) / 2;
                            if collision_len[k2][mid] <= *x {
                                a = mid + 1;
                            } else {
                                b = mid;
                            }
                        }
                        sum += a * x;
                        num_samples += a;
                    }

                    for x in &collision_len[k2] {
                        let (mut a, mut b) = (0, p);
                        while a < b {
                            let mid = (a + b) / 2;
                            if collision_len[k1][mid] < *x {
                                a = mid + 1;
                            } else {
                                b = mid;
                            }
                        }
                        sum += a * x;
                        num_samples += a;
                    }

                    assert_eq!(num_samples, p * p);

                    let expected = (sum as f64) / ((p * p) as f64);
                    println!(
                        "k1 = {:<4} k2 = {:<4} expected = {}",
                        k1, k2, expected
                    );
                }
            }
        }
    });
}
