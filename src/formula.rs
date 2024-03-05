use once_cell::sync::Lazy;

const K_MAX: u64 = 1 << 14;

struct Phi {
    phi: Vec<u64>,
    divisors: Vec<Vec<u64>>,
}

impl Phi {
    fn new(n: u64) -> Phi {
        let mut prime_factor = vec![0u64; n as usize + 1];

        for i in 2..=n {
            if prime_factor[i as usize] == 0 {
                prime_factor[i as usize] = i;
                let mut j = i * i;
                while j <= n {
                    prime_factor[j as usize] = i;
                    j += i;
                }
            }
        }

        let mut divisors: Vec<Vec<u64>> = vec![Vec::new(); n as usize + 1];

        for i in 1..=n {
            let mut j = i;
            while j <= n {
                divisors[j as usize].push(i);
                j += i;
            }
        }

        let mut phi: Vec<u64> = vec![0; n as usize + 1];
        phi[1] = 1;
        for i in 2..=n {
            let mut q = 1;
            let mut m = i;
            while m % prime_factor[i as usize] == 0 {
                q *= prime_factor[i as usize];
                m /= prime_factor[i as usize];
            }

            if m == 1 {
                phi[i as usize] = (prime_factor[i as usize] - 1)
                    * (q / prime_factor[i as usize]);
            } else {
                phi[i as usize] = phi[q as usize] * phi[(i / q) as usize];
            }
        }

        Phi { phi, divisors }
    }

    fn divisors(&self, n: u64) -> &Vec<u64> {
        &self.divisors[n as usize]
    }

    // The probability that gcd(n, x) = d for random 0 <= x < n. d must be a
    // divisor of n.
    fn gcd_probability(&self, n: u64, d: u64) -> f64 {
        self.phi[(n / d) as usize] as f64 / n as f64
    }
}

static PHI: Lazy<Phi> = Lazy::new(|| Phi::new(K_MAX));

fn independent_machines(k: &Vec<u64>, i: usize, s: f64) -> f64 {
    if i == k.len() {
        return 1.0 / s.sqrt();
    }

    // Silently assumes that the k-values are increasing. Search the last j
    // such that k[i] == k[j].
    let mut j = i;
    while j < k.len() && k[i] == k[j] {
        j += 1;
    }

    let inv_lg_k_i_squared =
        1.0 / (((k[i] << 1) as f64).log2() * ((k[i] << 1) as f64).log2());
    let mut expected: f64 = 0.0;
    for d in PHI.divisors(k[i]) {
        expected += independent_machines(
            k,
            j,
            s + (j as u64 * (2 * d - 1)) as f64 * inv_lg_k_i_squared,
        ) * PHI.gcd_probability(k[i], *d);
    }

    expected
}
// k is the array of k-values assigned to the machines.
pub fn expected_time(k: &Vec<u64>) -> f64 {
    independent_machines(k, 0, 0.0)
}
