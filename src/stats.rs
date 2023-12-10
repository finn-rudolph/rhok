use std::collections::VecDeque;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{miller_rabin::miller_rabin, montgomery::Montgomery, single::gcd};

const A: usize = 656159711;
const B: usize = 656159711;
const K1: usize = 3;
const K2: usize = 10;

fn f(x: usize, k: usize, mtg: &Montgomery) -> usize {
    mtg.strict(mtg.pow(x as u64, (k as u64) << 1) + 1) as usize
}

fn get_predecessors(p: usize, k: usize, mtg: &Montgomery) -> Vec<Vec<usize>> {
    let mut pre: Vec<Vec<usize>> = vec![Vec::new(); p];

    for i in 0..p {
        let next = f(i, k, mtg);
        pre[next].push(i);
    }

    pre
}

// iota: # predecessors that have in-degree > 1 (i.e lie in the interesting
// subgraph)
pub fn iota_stats() {
    println!(
        "{:<20}{:<16}{:<16}{:<20}{}",
        "", "k", "gcd(2k, p - 1)", "std dev(iota)", "avg of squares (iota)"
    );
    for p in A..=B {
        println!("p = {}", p,);
        let mtg = Montgomery::new(p as u64);

        for k in K1..=K2 {
            let pre = get_predecessors(p, k, &mtg);

            // mean is 1 since there is an equal number of in and out nodes.
            let mut var_iota: usize = 0;
            let mut quad_iota: usize = 0;
            let mut iota_histogram: [usize; 2 * K2 + 1] = [0; 2 * K2 + 1];
            let mut in_degree_histogram: [usize; 2 * K2 + 1] = [0; 2 * K2 + 1];

            for i in 0..p {
                // there's some weird in-degree 1 node, we wanna ignore this
                in_degree_histogram[pre[i].len()] += 1;
                if pre[i].len() > 0 {
                    let mut iota = 0;
                    for u in &pre[i] {
                        if pre[*u].len() > 1 {
                            iota += 1;
                        }
                    }
                    iota_histogram[iota] += 1;
                    quad_iota += iota * iota;
                    var_iota +=
                        ((iota as i64 - 1) * (iota as i64 - 1)) as usize;
                }
            }

            let g = gcd(p as u64 - 1, 2 * k as u64);
            println!(
                "{:<20}{:<16}{:<16}{:<20}{:<20} {:<?} {:?}",
                "",
                k,
                g,
                (var_iota as f64 / ((p as u64 / g) as f64)).sqrt(),
                (quad_iota as f64) / (p as u64 / g) as f64,
                iota_histogram,
                in_degree_histogram
            );

            // if g == 2 {
            //     if p & 7 == 1 {
            //         assert!(iota_histogram[0] == (p - 1) >> 3);
            //         assert!(iota_histogram[1] == (p + 3) >> 2);
            //         assert!(iota_histogram[2] == (p - 9) >> 3);
            //     } else if p & 7 == 3 {
            //         assert!(iota_histogram[0] == (p - 3) >> 3);
            //         assert!(iota_histogram[1] == (p + 1) >> 2);
            //         assert!(iota_histogram[2] == (p - 3) >> 3);
            //     } else if p & 7 == 5 {
            //         assert!(iota_histogram[0] == (p + 3) >> 3);
            //         assert!(iota_histogram[1] == (p - 1) >> 2);
            //         assert!(iota_histogram[2] == (p - 5) >> 3);
            //     } else if p & 7 == 7 {
            //         assert!(iota_histogram[0] == (p + 1) >> 3);
            //         assert!(iota_histogram[1] == (p - 3) >> 2);
            //         assert!(iota_histogram[2] == (p + 1) >> 3);
            //     }
            // }
        }
    }
}

// calculates average collision length (nu) for various k
pub fn mu_lambda_stats() {
    // println!(
    //     "{:<20}{:<16}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
    //     "",
    //     "k",
    //     "gcd(k, p - 1)",
    //     "nu mean",
    //     "nu std dev",
    //     "cycle mean",
    //     "cycle std dev",
    //     "tail mean",
    //     "tail std dev"
    // );

    let nu_means: Vec<[f64; K2 + 1]> = (A..=B)
        .into_par_iter()
        .map(|p| {
            let mut nu = [0.0; K2 + 1];

            if gcd(p as u64 - 1, 4) != 2
                || gcd(p as u64 - 1, 6) != 2
                || !miller_rabin(p as u64)
            {
                return nu;
            }

            // println!("p = {}", p,);

            let mtg = Montgomery::new(p as u64);
            let mut cycle_nodes: Vec<usize> = Vec::new();
            let mut collision_len: Vec<usize> = vec![0; p];
            let mut cycle_len: Vec<usize> = vec![0; p];

            for k in K1..=K2 {
                let pre = get_predecessors(p, k, &mtg);
                let mut visited: Vec<bool> = vec![false; p];

                let mut total_collision_len: usize = 0;
                let mut total_cycle_len: usize = 0;

                for i in 0..p {
                    if !visited[i] {
                        {
                            // mark cycle nodes
                            let (mut x, mut y) = (i, f(i, k, &mtg));
                            while x != y {
                                x = f(x, k, &mtg);
                                y = f(f(y, k, &mtg), k, &mtg);
                            }

                            cycle_nodes.clear();
                            while cycle_nodes.is_empty() || x != cycle_nodes[0]
                            {
                                cycle_nodes.push(x);
                                visited[x] = true;
                                x = f(x, k, &mtg);
                            }
                        }

                        for x in &cycle_nodes {
                            let mut q: VecDeque<(usize, usize)> =
                                VecDeque::new();
                            q.push_back((*x, 0));

                            while !q.is_empty() {
                                let (u, dis) = *q.front().unwrap();
                                q.pop_front();
                                visited[u] = true;
                                total_collision_len += dis + cycle_nodes.len();
                                collision_len[u] = dis + cycle_nodes.len();
                                cycle_len[u] = cycle_nodes.len();
                                total_cycle_len += cycle_nodes.len();

                                for v in &pre[u] {
                                    if !visited[*v] {
                                        q.push_back((*v, dis + 1));
                                    }
                                }
                            }
                        }
                    }
                }

                let mean_nu = (total_collision_len as f64) / (p as f64);
                nu[k] = mean_nu / (p as f64).sqrt();
                let mean_cycle_len = (total_cycle_len as f64) / (p as f64);
                let mean_tail_len = ((total_collision_len - total_cycle_len)
                    as f64)
                    / (p as f64);

                let mut std_dev_nu = 0.0;
                let mut std_dev_cycle_len = 0.0;
                let mut std_dev_tail_len = 0.0;
                for i in 0..p {
                    std_dev_nu += (collision_len[i] as f64 - mean_nu)
                        * (collision_len[i] as f64 - mean_nu);
                    std_dev_cycle_len += (cycle_len[i] as f64 - mean_cycle_len)
                        * (cycle_len[i] as f64 - mean_cycle_len);
                    std_dev_tail_len += ((collision_len[i] - cycle_len[i])
                        as f64
                        - mean_tail_len)
                        * ((collision_len[i] - cycle_len[i]) as f64
                            - mean_tail_len);
                }

                std_dev_nu = (std_dev_nu / p as f64).sqrt();
                std_dev_cycle_len = (std_dev_cycle_len / p as f64).sqrt();
                std_dev_tail_len = (std_dev_tail_len / p as f64).sqrt();

                // println!(
                //     "{:<20}{:<16}{:<16}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
                //     "",
                //     k,
                //     g,
                //     mean_nu,
                //     std_dev_nu,
                //     mean_cycle_len,
                //     std_dev_cycle_len,
                //     mean_tail_len,
                //     std_dev_tail_len
                // );
            }
            nu
        })
        .collect::<Vec<[f64; K2 + 1]>>()
        .into_iter()
        .filter(|x| x[2] != 0.0)
        .collect();

    let mut mean = [0.0; K2 + 1];
    let mut std_dev = [0.0; K2 + 1];
    for x in &nu_means {
        for k in K1..=K2 {
            mean[k] += x[k];
        }
    }
    for v in mean.iter_mut() {
        *v /= nu_means.len() as f64;
    }
    for x in &nu_means {
        for k in K1..=K2 {
            std_dev[k] += (x[k] - mean[k]) * (x[k] - mean[k]);
        }
    }
    for v in std_dev.iter_mut() {
        *v = (*v / nu_means.len() as f64).sqrt();
    }

    let mut corr = 0.0;
    let prob = 1.0 / nu_means.len() as f64;
    for x in &nu_means {
        corr += (x[2] - mean[2]) * (x[3] - mean[3]) * prob;
    }

    corr /= std_dev[2];
    corr /= std_dev[3];
    eprintln!("{}", corr);

    for x in &nu_means {
        println!("{}", x[2]);
    }

    for x in &nu_means {
        println!("{}", x[3]);
    }
}
