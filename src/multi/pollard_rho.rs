use rug::{rand::RandState, Complete, Integer};

fn pow_mod(mut x: Integer, mut y: u64, modulus: &Integer) -> Integer {
    while y & 1 == 0 {
        x = x.square() % modulus;
        y >>= 1;
    }
    if y == 1 {
        return x;
    }
    let mut result = x.clone();
    y ^= 1;

    loop {
        if y & 1 == 1 {
            result *= &x;
            result %= modulus;
        }
        y >>= 1;
        if y == 0 {
            break;
        }
        x = x.square() % modulus;
    }

    result
}

// TODO: set iter limit and use option

pub fn pollard_rho(
    n: &Integer,
    k: u64,
    rng: &mut RandState,
) -> Option<Integer> {
    const BATCH_SIZE: u64 = 1 << 10;
    const LENGTH_LIMIT: u64 = 1 << 18;

    let _k = k << 1;
    let c: Integer = n.random_below_ref(rng).into();

    let mut x: Integer = n.random_below_ref(rng).into();

    loop {
        let mut l = BATCH_SIZE;
        while l <= LENGTH_LIMIT {
            let (y, mut q) = (x.clone(), Integer::from(1));

            let mut i = 0;
            while i < l {
                for _ in 0..BATCH_SIZE {
                    x = pow_mod(x, _k, n) + &c;
                    q *= (&x - &y).complete();
                    q %= n;
                }

                let g = n.gcd_ref(&q).complete();
                if g != 1 && g != *n {
                    return Some(g);
                }

                i += BATCH_SIZE;
            }

            l <<= 1;
        }

        println!("retrying to factor {}", n);
    }
}
