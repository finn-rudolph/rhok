use rug::{
    rand::{RandGen, RandState},
    Complete, Integer,
};

// TODO: set iter limit and use option

pub fn pollard_rho(n: &Integer, k: u64, rng: &mut RandState) -> Integer {
    const BATCH_SIZE: u64 = 1 << 10;
    const LENGTH_LIMIT: u64 = 1 << 18;

    let _k = Integer::from(k << 1);

    loop {
        let mut x: Integer = n.random_below_ref(rng).into();

        let mut l = BATCH_SIZE;
        while l <= LENGTH_LIMIT {
            let (y, mut q) = (x.clone(), Integer::from(1));

            let mut i = 0;
            while i < l {
                for _ in 0..BATCH_SIZE {
                    x = x.pow_mod(&_k, n).unwrap();
                    q *= (&x - &y).complete();
                    q %= n;
                }

                let g = n.gcd_ref(&q).complete();
                if g != 1 && g != *n {
                    return g;
                }

                i += BATCH_SIZE;
            }

            l <<= 1;
        }
    }
}
