use crate::montgomery::Montgomery;

/// Determine if `p` is prime or not.
///
/// In the current version, the trait's method cannot be used in a const context,  
/// so this method can be used as a workaround only for `u64`.
///
/// This method will may be integrated into `IsPrime::is_prime` in a future update.
///
/// # Examples
/// ```rust
/// use primality_test::is_prime;
///
/// const PRIME: bool = is_prime(999999999999999989);
/// const NOT_PRIME: bool = is_prime(585226005592931977);
///
/// assert!(PRIME);
/// assert!(!NOT_PRIME);
/// ```
pub const fn is_prime(p: u64) -> bool {
    if p == 1 || p & 1 == 0 {
        return p == 2;
    }

    let mont = Montgomery::<u64>::new(p);

    let s = (p - 1).trailing_zeros();
    let t = (p - 1) >> s;

    let mut i = 0;
    const WIT: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    while i < WIT.len() {
        let a = WIT[i] % p;
        if a != 0 {
            let a = mont.convert(a);
            let mut at = mont.pow(a, t);

            // a^t = 1 (mod p) or a^t = -1 (mod p)
            let mut found = at == mont.r || at == p - mont.r;

            // found i satisfying a^((2^i)*t) = -1 (mod p)
            let mut j = 1;
            while j < s && !found {
                at = mont.multiply(at, at);
                found |= at == p - mont.r;
                j += 1;
            }

            if !found {
                return false;
            }
        }
        i += 1;
    }

    true
}

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
            assert!(is_prime(p))
        }
    }

    #[test]
    fn non_primality_test() {
        for &p in NO_PRIME {
            assert!(!is_prime(p));
        }

        for &p in CARMICHAEL {
            assert!(!is_prime(p));
        }
    }
}
