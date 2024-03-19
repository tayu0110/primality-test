#![allow(clippy::len_without_is_empty)]

use std::{
    array::IntoIter,
    iter::{from_fn, Enumerate, FilterMap},
};

/// Generate a table of prime numbers by linear sieving.
///
/// This structure supports table generation by compile-time calculations and can also determine prime numbers at compile time.  
///
/// # Limitations
/// The range generated is `0..LEN`, *NOT* `0..=LEN`. This is due to the limitations of the const generics.  
/// Therefore, even if memory and computation time were infinite, this structure could not be used to determine whether `usize::MAX` is prime or not.
///
/// Also, implementation of types other than usize is currently not provided.  
/// It may be possible to implement it through specialization or const_trait_impl stabilization, but no implementation method has been found so far.
///
/// # Examples
/// ```
/// use primality_test::LinearSieve;
///
/// const S: LinearSieve<15> = LinearSieve::new();
/// const PRIME: bool = S.is_prime(7);
/// const NOT_PRIME: bool = S.is_prime(10);
///
/// assert!(PRIME);
/// assert!(!NOT_PRIME);
/// ```
///
/// Note that `[T; N]` is used instead of `Vec<T>` to support compile-time calculations.  
/// This may result in increased stack memory usage for runtime table generation.
/// ```no_run
/// use primality_test::LinearSieve;
///
/// // Data size is too big...
/// // At this line, stack overflow may occur.
/// let _s = LinearSieve::<10000000000000>::new();
/// ```
#[derive(Debug, Clone, Copy)]
pub struct LinearSieve<const LEN: usize> {
    arr: [usize; LEN],
}

impl<const LEN: usize> LinearSieve<LEN> {
    /// Constructor.  
    /// This method supports the compile-time generation.
    pub const fn new() -> Self {
        let mut primes = [0; LEN];
        let mut pc = 0;
        let mut arr = [usize::MAX; LEN];
        let mut i = 2;
        while i < LEN {
            if arr[i] == usize::MAX {
                primes[pc] = i;
                pc += 1;
                arr[i] = i;
            }
            let mut j = 0;
            while j < pc && primes[j] <= arr[i] && primes[j].saturating_mul(i) < LEN {
                arr[primes[j] * i] = primes[j];
                j += 1;
            }
            i += 1;
        }

        Self { arr }
    }

    /// Return the length of the internal array.  
    ///
    /// The value returned by this method is always equal to `LEN`.
    pub const fn len(&self) -> usize {
        self.arr.len()
    }

    /// Check if `value` is prime or not.
    ///
    /// # Panics
    /// This method panics when `value >= LEN` is `true`.
    ///
    /// # Examples
    /// ```
    /// use primality_test::LinearSieve;
    ///
    /// const S: LinearSieve<1000> = LinearSieve::new();
    /// assert!(S.is_prime(2));
    /// assert!(S.is_prime(997));
    ///
    /// assert!(!S.is_prime(0));
    /// assert!(!S.is_prime(1));
    /// assert!(!S.is_prime(561));
    /// // S.is_prime(1000)     // Panic at this line
    /// ```
    pub const fn is_prime(&self, value: usize) -> bool {
        self.arr[value] == value
    }

    /// Factorize `value` into prime factors.
    ///
    /// `LinearSieve` keeps internally the least prime factor of each value.  
    /// By using this, the method is guaranteed to return the prime factors of `value` in ascending order.
    ///
    /// `1` is not contained in the returned factors.  
    /// If `value` is less then `2`, the length of returned iterator is always `0`.
    ///
    /// # Panics
    /// This method panics when `value` >= `LEN` is `true`.
    ///
    /// # Examples
    /// ```
    /// use primality_test::LinearSieve;
    ///
    /// const S: LinearSieve<100> = LinearSieve::new();
    /// let factors_of_120 = S.factors(60).collect::<Vec<_>>();
    /// assert_eq!(factors_of_120, vec![2, 2, 3, 5]);
    ///
    /// let factors_of_11 = S.factors(11).collect::<Vec<_>>();
    /// assert_eq!(factors_of_11, vec![11]);
    ///
    /// assert!(S.factors(1).next().is_none());
    /// ```
    pub fn factors(&self, mut value: usize) -> impl Iterator<Item = usize> + '_ {
        from_fn(move || {
            (value > 1).then(|| {
                let res = self.arr[value];
                value /= res;
                res
            })
        })
    }
}

impl<const LEN: usize> IntoIterator for LinearSieve<LEN> {
    type Item = usize;
    type IntoIter = FilterMap<Enumerate<IntoIter<usize, LEN>>, fn((usize, usize)) -> Option<usize>>;
    /// Generate a stream of prime numbers.
    fn into_iter(self) -> Self::IntoIter {
        self.arr
            .into_iter()
            .enumerate()
            .filter_map(|(i, p)| (i == p).then_some(p))
    }
}

impl<const LEN: usize> Default for LinearSieve<LEN> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Mul;

    use crate::*;

    const LEN: usize = 10000;
    const S: LinearSieve<LEN> = LinearSieve::new();

    #[test]
    fn prime_check_test() {
        for i in 0..LEN {
            assert_eq!(i.is_prime(), S.is_prime(i));
        }
    }

    #[test]
    fn sieve_iterator_test() {
        assert!(S.into_iter().all(|p| p.is_prime()))
    }

    #[test]
    fn factors_test() {
        for i in 0..LEN {
            let f = S.factors(i).collect::<Vec<_>>();
            assert!(f.windows(2).all(|v| v[0] <= v[1]));
            if i > 1 {
                assert_eq!(f.into_iter().reduce(Mul::mul).unwrap(), i);
            } else {
                assert!(f.is_empty());
            }
        }
    }
}
