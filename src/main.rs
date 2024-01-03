use rhok::{miller_rabin::miller_rabin, montgomery::Montgomery};

fn main() {
    rhok::single::bench_single_rho();
    // for p in 1 << 22..1 << 23 {
    //     if (p - 1) & 3 != 0 || !miller_rabin(p) {
    //         continue;
    //     }

    //     let mtg = Montgomery::new(p);
    //     let mut cnt = 0;
    //     let mut w = 0;
    //     let neg_one = mtg.to_montgomery_space(p - 1);
    //     let mut z = neg_one;
    //     for _ in 1..p / 2 {
    //         if mtg.strict(mtg.pow(w, (p - 1) >> 2)) == mtg.one()
    //             && mtg.strict(mtg.pow(z, (p - 1) >> 2)) == mtg.one()
    //         {
    //             cnt += 1;
    //         }
    //         w = mtg.add(w, mtg.one());
    //         z = mtg.add(z, neg_one);
    //     }
    //     println!("c = {}, (p - 1) / 32 = {}", cnt, (p - 1) >> 5);
    // }
}
