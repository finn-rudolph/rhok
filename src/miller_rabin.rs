use crate::montgomery::Montgomery;
// Deterministic variant of Miller-Rabin that works up to 2^64.
pub const fn miller_rabin(n: u64) -> bool {
    if n == 2 {
        return true;
    }
    if n & 1 == 0 {
        return false;
    }

    const MILLER_RABIN_BASES: [u64; 7] =
        [2, 325, 9375, 28178, 450775, 9780504, 1795265022];

    let trailing_zeros = (n - 1).trailing_zeros();
    let u = (n - 1) >> trailing_zeros;

    let mtg = Montgomery::new(n);
    let one_mtg = mtg.strict(mtg.to_montgomery_space(1));
    let neg_one_mtg = mtg.strict(mtg.to_montgomery_space(n - 1));

    let mut i = 0;
    while i < MILLER_RABIN_BASES.len() {
        let mut a = MILLER_RABIN_BASES[i];
        a %= n;
        if a == 0 {
            i += 1;
            continue;
        }
        a = mtg.to_montgomery_space(a);

        let mut x = mtg.strict(mtg.pow(a, u));
        let mut j = 0;
        while j < trailing_zeros {
            let y = mtg.strict(mtg.mul(x, x));
            if y == one_mtg && x != one_mtg && x != neg_one_mtg {
                return false;
            }
            x = y;
            j += 1;
        }
        if x != one_mtg {
            return false;
        }
        i += 1;
    }

    true
}

#[cfg(test)]
mod test {
    use rand_xoshiro::{
        rand_core::{RngCore, SeedableRng},
        Xoshiro256PlusPlus,
    };
    use rug::{integer::IsPrime, Integer};

    use super::*;

    #[test]
    fn test_miller_rabin() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(71 * 42);

        for _ in 0..400000 {
            let n = rng.next_u64() >> 2;
            assert_eq!(
                miller_rabin(n),
                Integer::from(n).is_probably_prime(64) != IsPrime::No
            );
        }

        for n in (2..400000)
            .chain([2, 3, 5, 13, 19, 73, 193, 407521, 299210837].into_iter())
        {
            println!("{}", n);
            assert_eq!(
                miller_rabin(n),
                Integer::from(n).is_probably_prime(64) != IsPrime::No
            );
        }
    }
}
