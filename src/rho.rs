use rug::{rand::RandState, Integer};

const BATCH_SIZE: u64 = 1 << 9;
const LENGTH_LIMIT: u64 = 1 << 17;

pub fn pollard_rho(n: &Integer, k: u32) -> Integer {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(42));

    let (mut x, mut y, mut q, mut d) = (
        Integer::new(),
        Integer::new(),
        Integer::new(),
        Integer::new(),
    );
    let _k = Integer::from(k);

    loop {
        x = n.random_below_ref(&mut rng).into();

        let mut l = BATCH_SIZE;
        while l <= LENGTH_LIMIT {
            y = x.clone();
            q = 1.into();

            let mut i = 0;
            while i < l {
                for _ in 0..BATCH_SIZE {
                    x = x.pow_mod(&_k, &n).unwrap() + 1;
                    d = (&x - &y).into();
                    q *= d;
                }

                d = q.gcd_ref(n).into();
                if d != 1 && d != *n {
                    return d;
                }

                i += BATCH_SIZE;
            }

            l <<= 1;
        }
    }
}
