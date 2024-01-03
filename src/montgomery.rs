pub(crate) struct Montgomery<T> {
    pub(crate) modulo: T,
    pub(crate) modulo_inv: T,
    pub(crate) r: T,
    pub(crate) r2: T,
}

macro_rules! impl_primitive_montgomery {
    ( $t:ty, $expand:ty ) => {
        impl Montgomery<$t> {
            // t <- MR(T) = floor(T/R) - floor((TN' mod R)*N/R)
            //  if t < 0 then return t + N else return t
            //      T := a (0 <= T < NR)
            //      N := MOD
            //      N':= MOD_INV    NN' = 1 (mod R)
            //      R := R
            #[allow(unused)]
            pub const fn reduce(&self, val: $t) -> $t {
                let (t, f) = (((val.wrapping_mul(self.modulo_inv) as $expand)
                    .wrapping_mul(self.modulo as $expand)
                    >> <$t>::BITS) as $t)
                    .overflowing_neg();
                t.wrapping_add(self.modulo * f as $t)
            }

            pub const fn multiply(&self, lhs: $t, rhs: $t) -> $t {
                let a = lhs as $expand * rhs as $expand;
                let (t, f) = ((a >> <$t>::BITS) as $t).overflowing_sub(
                    (((a as $t).wrapping_mul(self.modulo_inv) as $expand)
                        .wrapping_mul(self.modulo as $expand)
                        >> <$t>::BITS) as $t,
                );
                t.wrapping_add(self.modulo * f as $t)
            }

            pub const fn convert(&self, val: $t) -> $t {
                self.multiply(val, self.r2)
            }

            pub const fn pow(&self, val: $t, mut exp: $t) -> $t {
                let (mut res, mut val) = (self.r, val);
                while exp > 0 {
                    if exp & 1 != 0 {
                        res = self.multiply(res, val);
                    }
                    val = self.multiply(val, val);
                    exp >>= 1;
                }
                res
            }

            pub const fn new(modulo: $t) -> Self {
                let r = (((1 as $expand) << <$t>::BITS) % modulo as $expand) as $t;
                let r2 = ((modulo as $expand).wrapping_neg() % modulo as $expand) as $t;
                let modulo_inv = {
                    let mut inv = modulo;
                    while modulo.wrapping_mul(inv) != 1 {
                        inv = inv.wrapping_mul((2 as $t).wrapping_sub(modulo.wrapping_mul(inv)));
                    }
                    inv
                };
                Self {
                    modulo,
                    modulo_inv,
                    r,
                    r2,
                }
            }
        }
    };
}

impl_primitive_montgomery!(u8, u16);
impl_primitive_montgomery!(u16, u32);
impl_primitive_montgomery!(u32, u64);
impl_primitive_montgomery!(u64, u128);
