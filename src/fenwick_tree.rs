use std::ops::{AddAssign, Sub};

pub struct FenwickTree<T>(Vec<T>);

impl<'a, T> FenwickTree<T>
where
    T: Default + Copy + AddAssign<T> + Sub<Output = T>,
{
    pub fn new(n: usize) -> FenwickTree<T> {
        FenwickTree(vec![Default::default(); n])
    }

    pub fn update(&mut self, i: usize, delta: T) {
        let mut j = i as isize + 1;
        while j != self.0.len() as isize {
            self.0[j as usize - 1] += delta;
            j += j & -j;
        }
    }

    pub fn prefix_sum(&self, i: usize) -> T {
        let mut x: T = Default::default();
        let mut j = i as isize + 1;
        while j != 0 {
            x += self.0[j as usize - 1];
            j -= j & -j;
        }
        x
    }

    pub fn range_sum(&self, i: usize, j: usize) -> T {
        if i == 0 {
            self.prefix_sum(j)
        } else {
            self.prefix_sum(j) - self.prefix_sum(i - 1)
        }
    }
}
