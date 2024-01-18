use core::cmp;

use rand::RngCore;

use crate::montgomery::Montgomery;

pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }

    let mut a_trz = a.trailing_zeros();
    let b_trz = b.trailing_zeros();
    let s = cmp::min(a_trz, b_trz);
    b >>= s;

    while a != 0 {
        a >>= a_trz;
        let d = b.wrapping_sub(a) as i64;
        a_trz = d.trailing_zeros();
        b = cmp::min(a, b);
        a = (if d < 0 { -d } else { d }) as u64;
    }

    b << s
}

// Pollard's Rho algorithm with Brent's optimization, taken from Sergey Slotin's
// book "Algorithms for Modern Hardware". (the gcd above is also from there)
pub fn pollard_rho(n: u64, k: u64, rng: &mut dyn RngCore) -> u64 {
    const BATCH_SIZE: u64 = 1 << 10;
    const LENGTH_LIMIT: u64 = 1 << 18;

    let mtg = Montgomery::new(n);
    let one = mtg.to_montgomery_space(1);
    let _k = k << 1;

    let mut x = rng.next_u64() % n;

    let mut l = BATCH_SIZE;
    while l <= LENGTH_LIMIT {
        let (y, mut q) = (x, one);

        let mut i = 0;
        while i < l {
            for _ in 0..BATCH_SIZE {
                x = mtg.add(mtg.pow(x, _k), one);
                let mut d = x.wrapping_sub(y);
                d = d.wrapping_add(if (d as i64) < 0 { n } else { 0 });
                q = mtg.mul(q, d);
            }

            let g = gcd(q, n);
            if g != 1 && g != n {
                return g;
            }

            i += BATCH_SIZE;
        }

        l <<= 1;
    }

    return u64::MAX;
}
