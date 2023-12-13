pub struct Montgomery {
    pub n: u64, // the modulus
    two_n: u64, // 2 * n
    n_inv: u64, // n^-1 mod 2^64
    one: u64,
    two_to_128_mod_n: u64, // 2^128 mod n
}

impl Montgomery {
    pub const fn new(n: u64) -> Montgomery {
        let mut n_inv: u64 = 1;
        let mut i = 0;
        while i < 6 {
            n_inv =
                n_inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(n_inv)));
            i += 1;
        }
        let mut mtg = Montgomery {
            n,
            two_n: n << 1,
            n_inv,
            one: 0,
            two_to_128_mod_n: {
                let two_to_64_mod_n = (1u128 << 64) % n as u128;
                ((two_to_64_mod_n * two_to_64_mod_n) % n as u128) as u64
            },
        };
        mtg.one = mtg.to_montgomery_space(1);
        mtg
    }

    pub const fn n(&self) -> u64 {
        self.n
    }

    #[inline(always)]
    pub const fn to_montgomery_space(&self, x: u64) -> u64 {
        self.mul(x, self.two_to_128_mod_n)
    }

    // Divides t by 2^64, like in the regular reduction, but with minor
    // optimizations possible since the input is only 64-bit.
    #[inline(always)]
    pub const fn out_of_montgomery_space(&self, t: u64) -> u64 {
        let m = t.wrapping_mul(self.n_inv);
        let mn_high = ((m as u128 * self.n as u128) >> 64) as u64;
        // we know mn_high <= n and t_high is zero here
        self.n - mn_high
    }

    #[inline(always)]
    pub const fn add(&self, x: u64, y: u64) -> u64 {
        let z = x.wrapping_add(y);
        z - if z >= self.two_n { self.two_n } else { 0 }
    }

    // Computes the product of x and y in Montgomery Space, reducing the result
    // to 0..2n - 1.
    #[inline(always)]
    pub const fn mul(&self, x: u64, y: u64) -> u64 {
        let t = x as u128 * y as u128;
        let (t_low, t_high) = (t as u64, (t >> 64) as u64);

        // m = t * n^-1 mod 2^64
        let m = t_low.wrapping_mul(self.n_inv);
        // We don't need to consider the lower bits of mn, since they will
        // cancel anyway in the following subtraction (since t - mn is
        // divisible by 2^64).
        let mn_high = ((m as u128 * self.n as u128) >> 64) as u64;
        t_high + self.n - mn_high
    }

    // x^y, reduced to 0..2n - 1. Of course, only x should be in Montgomery
    // Space!
    #[inline]
    pub const fn pow(&self, mut x: u64, mut y: u64) -> u64 {
        let mut result: u64 = x;
        while y & 1 == 0 {
            result = self.mul(result, result);
            y >>= 1;
        }
        x = result;
        y ^= 1;

        loop {
            if y & 1 == 1 {
                result = self.mul(result, x);
            }
            y >>= 1;
            if y == 0 {
                break;
            }
            x = self.mul(x, x);
        }

        result
    }

    #[inline(always)]
    pub const fn strict(&self, x: u64) -> u64 {
        x - if x >= self.n { self.n } else { 0 }
    }

    #[inline(always)]
    pub const fn one(&self) -> u64 {
        self.one
    }
}
