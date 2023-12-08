use crate::{miller_rabin::miller_rabin, montgomery::Montgomery, single::gcd};

const A: usize = 1 << 17;
const B: usize = (1 << 17) + (1 << 15);
const K1: usize = 2;
const K2: usize = 12;

pub fn stats() {
    println!(
        "{:<20}{:<16}{:<16}{:<20}{}",
        "", "k", "gcd(k, p - 1)", "std dev(iota)", "avg of squares (iota)"
    );
    for p in A..=B {
        if !miller_rabin(p as u64) {
            continue;
        }

        println!("p = {}", p,);
        let mtg = Montgomery::new(p as u64);
        let mut in_degree: Vec<usize> = vec![0; p];
        let mut pre: Vec<[usize; K2 + 1]> = vec![[0; K2 + 1]; p];

        for k in (K1..=K2).step_by(2) {
            in_degree.fill(0);
            for i in 0..p {
                let next = mtg.strict(mtg.pow(i as u64, k as u64) + 1) as usize;
                pre[next][in_degree[next]] = i;
                in_degree[next] += 1;
            }

            // mean is 1 since there is an equal number of in and out nodes.
            let mut var_iota: usize = 0;
            let mut quad_iota: usize = 0;
            let mut histogram: [usize; K2 + 1] = [0; K2 + 1];
            for i in 0..p {
                // there's some weird in-degree 1 node, we wanna ignore this
                if in_degree[i] > 1 {
                    let mut iota = 0;
                    for j in 0..in_degree[i] {
                        if in_degree[pre[i][j]] > 1 {
                            iota += 1;
                        }
                    }
                    histogram[iota] += 1;
                    quad_iota += iota * iota;
                    var_iota +=
                        ((iota as i64 - 1) * (iota as i64 - 1)) as usize;
                }
            }

            let g = gcd(p as u64 - 1, k as u64);
            println!(
                "{:<20}{:<16}{:<16}{:<20}{:<20} {:<?}",
                "",
                k,
                g,
                (var_iota as f64 / ((p as u64 / g) as f64)).sqrt(),
                (quad_iota as f64) / (p as u64 / g) as f64,
                histogram
            );
        }
    }
}
