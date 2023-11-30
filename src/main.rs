use std::env;

use rug::{Complete, Integer};

fn main() {
    let args: Vec<String> = env::args().collect();

    let n = Integer::parse(&args[1])
        .expect("failed to parse the number to be factored")
        .complete();

    println!("=> {}", rhok::pollard_rho(&n, 2));
}
