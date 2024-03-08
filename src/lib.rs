//! Provide a method to determine whether an unsigned integer is a prime number through the `IsPrime` Trait.
//!
//! In the current implementation, this crate uses the Miller-Rabin primality test.  
//! The Miller-Rabin primality test is known to have witnesses that can conclusively determine unsigned integers of at most 64 bits.  
//! In this crate, the following information is used to select the witnesses.  
//!
//! [Deterministic variants of the Miller-Rabin primality test](https://miller-rabin.appspot.com/)

mod montgomery;
mod sieve;

use montgomery::Montgomery;
pub use sieve::LinearSieve;

const SMALL_PRIMES_MEMO: LinearSieve<255> = LinearSieve::new();

pub trait IsPrime {
    fn is_prime(&self) -> bool;
}

macro_rules! impl_is_prime {
    ( $t:ty, $witness_ty:ty, [ $( $witness:expr ),* ] ) => {
        impl IsPrime for $t {
            /// Determine an unsigned integer is prime number or not.
            ///
            /// # Examples
            /// ```rust
            /// use primality_test::IsPrime;
            ///
            /// assert!(998244353u32.is_prime());
            /// assert!(!561u16.is_prime());
            ///
            /// let primes = (1..20u16).filter(IsPrime::is_prime).collect::<Vec<_>>();
            /// assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
            /// ```
            fn is_prime(&self) -> bool {
                let p = *self as $witness_ty;

                if p == 1 || p & 1 == 0 {
                    return p == 2;
                }

                if <$witness_ty>::BITS == 64 && (p as u64) < 1795265022 {
                    return (p as u32).is_prime();
                }

                if p < SMALL_PRIMES_MEMO.len() as $witness_ty {
                    return SMALL_PRIMES_MEMO.is_prime(p as usize);
                }

                let mont = Montgomery::<$witness_ty>::new(p);

                let s = (p - 1).trailing_zeros();
                let t = (p - 1) >> s;

                [$( $witness ),*]
                    .iter()
                    // These two lines are necessary if any of the witnesses may be greater than or equal to `p`.
                    // The maximum witness for upto 32-bit integers is 61, and when `p` <= 61, this line is unreachable because they dictionary SMALL_PRIMES_MEMO.
                    // The maximum witness for 64-bit integers is 1795265022, but since numbers less than this can be expressed in 32 bits, they are cast to u32 and judged again, so this point is unreachable.
                    // .map(|&a| a % p)
                    // .filter(|&a| a != 0)
                    .all(|&a| {
                        let a = mont.convert(a);
                        let at = mont.pow(a, t);
                        // a^t = 1 (mod p) or a^t = -1 (mod p)
                        if at == mont.r || at == p - mont.r {
                            return true;
                        }

                        // found i satisfying a^((2^i)*t) = -1 (mod p)
                        (1..s)
                            .scan(at, |at, _| {
                                *at = mont.multiply(*at, *at);
                                Some(*at)
                            })
                            .any(|at| at == p - mont.r)
                    })
            }
        }
    };
}

impl_is_prime!(u8, u8, [2, 7, 61]);
impl_is_prime!(u16, u16, [2, 7, 61]);
impl_is_prime!(u32, u32, [2, 7, 61]);
impl_is_prime!(u64, u64, [2, 325, 9375, 28178, 450775, 9780504, 1795265022]);
#[cfg(target_pointer_width = "64")]
impl_is_prime!(
    usize,
    u64,
    [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
);
#[cfg(target_pointer_width = "32")]
impl_is_prime!(usize, u32, [2, 7, 61]);

#[cfg(test)]
mod tests {
    use super::*;

    #[rustfmt::skip]
    const PRIME: &[u64] = &[
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 73, 79, 83, 89, 97,
        998244353, 1000000007, 67280421310721, 999999999999999989
    ];

    const NO_PRIME: &[u64] = &[
        1,
        57,
        5329,
        49141,
        4759123141,
        1122004669633,
        21652684502221,
        31858317218647,
        47636622961201,
        55245642489451,
        3071837692357849,
        3770579582154547,
        7999252175582851,
        585226005592931977,
    ];

    #[rustfmt::skip]
    const CARMICHAEL: &[u64] = &[
        561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041, 46657, 52633, 62745,
        63973, 75361, 101101, 115921, 126217, 162401, 172081, 188461, 252601, 278545, 294409,
        314821, 334153, 340561, 399001, 410041, 449065, 488881, 512461,
    ];

    #[test]
    fn primality_test() {
        for &p in PRIME {
            if p <= u8::MAX as u64 {
                assert!((p as u8).is_prime());
            }
            if p <= u16::MAX as u64 {
                assert!((p as u16).is_prime());
            }
            if p <= u32::MAX as u64 {
                assert!((p as u32).is_prime());
            }
            assert!(p.is_prime());
        }
    }

    #[test]
    fn non_primality_test() {
        for &p in NO_PRIME {
            if p <= u8::MAX as u64 {
                assert!(!(p as u8).is_prime());
            }
            if p <= u16::MAX as u64 {
                assert!(!(p as u16).is_prime());
            }
            if p <= u32::MAX as u64 {
                assert!(!(p as u32).is_prime());
            }
            assert!(!p.is_prime());
        }

        for &p in CARMICHAEL {
            if p <= u8::MAX as u64 {
                assert!(!(p as u8).is_prime());
            }
            if p <= u16::MAX as u64 {
                assert!(!(p as u16).is_prime());
            }
            if p <= u32::MAX as u64 {
                assert!(!(p as u32).is_prime());
            }
            assert!(!p.is_prime());
        }
    }
}
