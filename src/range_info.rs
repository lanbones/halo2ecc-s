use halo2_proofs::arithmetic::BaseExt;
use halo2_proofs::arithmetic::FieldExt;
use num_bigint::BigUint;
use num_integer::Integer;
use std::marker::PhantomData;

use crate::utils::{bn_to_field, field_to_bn};

pub const MAX_CHUNKS: u64 = 5;
pub const MAX_BITS: u64 = 18;
pub const COMMON_RANGE_BITS: u64 = MAX_BITS as u64;

// 3 LIMBS
pub const W_CEIL_LEADING_CHUNKS: u64 = 5;
pub const N_FLOOR_LEADING_CHUNKS: u64 = 5;
pub const D_LEADING_CHUNKS: u64 = 5;
pub const LIMBS: usize = 3;

// 4 LIMBS
/*
pub const W_CEIL_LEADING_CHUNKS: u64 = 3;
pub const N_FLOOR_LEADING_CHUNKS: u64 = 3;
pub const D_LEADING_CHUNKS: u64 = 4;
pub const MAX_CHUNKS: u64 = 4;
pub const LIMBS: u64 = 4;
*/

pub const LIMB_BITS: u64 = MAX_CHUNKS as u64 * COMMON_RANGE_BITS;

pub const OVERFLOW_BITS: u64 = 6;
pub const OVERFLOW_LIMIT: u64 = 1 << OVERFLOW_BITS;
pub const OVERFLOW_THRESHOLD: u64 = 1 << (OVERFLOW_BITS - 2);

#[derive(Debug, Clone)]
pub struct RangeInfo<N: FieldExt> {
    pub w_ceil_bits: u64,
    pub d_leading_bits: u64,
    pub w_ceil_leading_bits: u64,
    pub n_floor_leading_bits: u64,

    pub w_ceil: BigUint,
    pub n_modulus: BigUint,
    pub neg_w_modulus: BigUint,
    pub w_modulus: BigUint,
    pub common_range_mask: BigUint,
    pub limb_mask: BigUint,
    pub limb_modulus: BigUint,
    pub integer_modulus: BigUint,
    pub max_d: BigUint,

    pub w_modulus_limbs_le_bn: [BigUint; LIMBS],
    pub w_modulus_limbs_le: [N; LIMBS],
    pub neg_w_modulus_limbs_le: [N; LIMBS],
    pub limb_coeffs: [N; LIMBS],
    pub limb_modulus_n: N,

    pub w_native: N,

    pub _phantom: PhantomData<N>,
}

impl<N: FieldExt> RangeInfo<N> {
    pub fn pre_check(&self) {
        // is_pure_w_modulus():
        // lcm(limb, native) >= w_ceil
        let limb_modulus = &self.limb_modulus;
        let lcm = self.n_modulus.lcm(&limb_modulus);
        assert!(lcm >= self.w_ceil);

        // reduce():
        // lcm(limb, native) >= w_ceil * OVERFLOW_LIMIT
        let limb_modulus = &self.limb_modulus;
        let lcm = self.n_modulus.lcm(&limb_modulus);
        assert!(lcm >= &self.w_ceil * OVERFLOW_LIMIT);
        // to ensure that d can be assigned by assign_common.
        assert!(COMMON_RANGE_BITS > OVERFLOW_BITS);
        // to ensure that v can be assigned by assign_nonleading_limb
        assert!(
            &(BigUint::from(1u64) << COMMON_RANGE_BITS + 1 + OVERFLOW_BITS) < &self.limb_modulus
        );
        // to ensure that v * limb_modulus would not overflow
        assert!(&(&self.limb_modulus * &self.limb_modulus) < &self.n_modulus);

        // mul():
        // lcm(integer_modulus, native) >= w_ceil * w_ceil * OVERFLOW_LIMIT * OVERFLOW_LIMIT
        let lcm = self.integer_modulus.lcm(&self.n_modulus);
        let max_a = &self.w_modulus * OVERFLOW_LIMIT;
        let max_b = &self.w_modulus * OVERFLOW_LIMIT;
        let max_l = max_a * max_b;
        let max_d = &self.max_d;
        let max_w = &self.w_modulus;
        let max_rem = &self.w_ceil;
        let max_r = max_d * max_w + max_rem;
        assert!(max_l <= lcm);
        assert!(max_r <= lcm);
        assert!(max_l <= max_r);
        // (u, v)
        let sum_limb_max = LIMBS as u64
            * (OVERFLOW_LIMIT * OVERFLOW_LIMIT + 1)
            * &self.limb_modulus
            * &self.limb_modulus;
        assert!(sum_limb_max < self.n_modulus);
        // Ensure we can split v into (common, nonleading)
        // The last `1` reserve for the borrow and carry.
        // The real reserved value is (1 + common) / limb_modulus.
        assert!(1 << COMMON_RANGE_BITS > LIMBS as u64 * (OVERFLOW_LIMIT * OVERFLOW_LIMIT + 1 + 1));
        let common_modulus = BigUint::from(1u64) << COMMON_RANGE_BITS;
        assert!(&common_modulus + 1u64 < self.limb_modulus);

        // Ensure nonoverflow for one limb sum check.
        assert!(&self.limb_modulus * &self.limb_modulus * common_modulus < self.n_modulus);

        // algorithm limitation
        assert!(LIMBS >= 3);
    }

    pub fn d_bits<W: BaseExt>() -> u64 {
        // a * b = w * d + rem

        // -- a <= w_ceil_bits + overflow_bits /\ b <= w_ceil_bits + overflow_bits
        // -> a * b <= w_ceil_bits * 2 + overflow_bits * 2 - 1

        // -- w * d + rem >= a * b
        // <- d_bits + w_ceil_bits - 2 = w_ceil_bits * 2 + overflow_bits * 2 - 1

        let w_max = field_to_bn(&-W::one());
        let w_ceil_bits = w_max.bits();
        let d_bits = w_ceil_bits + OVERFLOW_BITS * 2 + 1;
        assert!(d_bits + w_ceil_bits - 2 >= w_ceil_bits * 2 + OVERFLOW_BITS * 2 - 1);
        d_bits
    }

    pub fn new<W: BaseExt>() -> Self {
        let w_max = field_to_bn(&-W::one());
        let w_ceil_bits = w_max.bits();
        assert!(BigUint::from(1u64) << w_ceil_bits > w_max);
        assert!(BigUint::from(1u64) << (w_ceil_bits - 1) < w_max);
        let w_ceil_leading_bits = w_ceil_bits % COMMON_RANGE_BITS as u64;
        // limited by current implementation
        assert_eq!(
            w_ceil_bits / COMMON_RANGE_BITS % MAX_CHUNKS,
            W_CEIL_LEADING_CHUNKS as u64 - 1
        );
        assert!(w_ceil_leading_bits != 0);

        let n_max = field_to_bn(&-N::one());
        let n_floor_bits = n_max.bits() - 1;
        assert!(BigUint::from(1u64) << n_floor_bits < n_max);
        assert!(BigUint::from(1u64) << (n_floor_bits + 1) >= n_max);
        let n_floor_leading_bits = n_floor_bits % COMMON_RANGE_BITS;
        // limited by current implementation
        assert_eq!(
            n_floor_bits / COMMON_RANGE_BITS % MAX_CHUNKS,
            N_FLOOR_LEADING_CHUNKS - 1
        );
        assert!(n_floor_leading_bits != 0);

        let d_range_bits = Self::d_bits::<W>();
        let d_leading_bits = d_range_bits % COMMON_RANGE_BITS;
        // limited by current implementation
        assert_eq!(
            d_range_bits / COMMON_RANGE_BITS % MAX_CHUNKS,
            D_LEADING_CHUNKS - 1
        );
        assert!(d_leading_bits != 0);

        let max_d = BigUint::from(1u64) << d_range_bits;
        let limb_mask = (BigUint::from(1u64) << LIMB_BITS) - 1u64;
        let n_modulus = &n_max + 1u64;
        let w_modulus = &w_max + 1u64;
        let w_native = &w_modulus % &n_modulus;
        let integer_modulus = BigUint::from(1u64) << (LIMBS as u64 * LIMB_BITS);
        let neg_w_modulus = &integer_modulus - &w_modulus;

        let w_modulus_limbs_le_bn: [BigUint; LIMBS] = (0..LIMBS as u64)
            .into_iter()
            .map(|i| ((&w_modulus >> (i * LIMB_BITS)) & &limb_mask))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let neg_w_modulus_limbs_le: [N; LIMBS] = (0..LIMBS as u64)
            .into_iter()
            .map(|i| (bn_to_field::<N>(&((&neg_w_modulus >> (i * LIMB_BITS)) & &limb_mask))))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let w_modulus_limbs_le = w_modulus_limbs_le_bn
            .iter()
            .map(|x| bn_to_field(x))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();

        let limb_modulus = BigUint::from(1u64) << LIMB_BITS;
        let limb_modulus_n = bn_to_field(&limb_modulus);

        let res = Self {
            w_ceil_bits,
            d_leading_bits,
            w_ceil_leading_bits,
            n_floor_leading_bits,

            max_d,
            w_ceil: BigUint::from(1u64) << w_ceil_bits,
            n_modulus,
            w_modulus,
            neg_w_modulus,
            common_range_mask: BigUint::from((1u64 << COMMON_RANGE_BITS) - 1),
            limb_mask,
            limb_modulus,
            integer_modulus,
            limb_coeffs: (0..LIMBS as u64)
                .into_iter()
                .map(|i| bn_to_field(&(BigUint::from(1u64) << (i * LIMB_BITS))))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),

            limb_modulus_n,
            neg_w_modulus_limbs_le,
            w_modulus_limbs_le,
            w_modulus_limbs_le_bn,

            w_native: bn_to_field::<N>(&w_native),

            _phantom: PhantomData,
        };

        res.pre_check();
        res
    }

    pub fn bn_to_limb_le_n(&self, w: &BigUint) -> [N; LIMBS] {
        (0..LIMBS as u64)
            .into_iter()
            .map(|i| bn_to_field(&((w >> (i * LIMB_BITS)) & &self.limb_mask)))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }

    pub fn bn_to_limb_le(&self, w: &BigUint) -> [BigUint; LIMBS] {
        (0..LIMBS as u64)
            .into_iter()
            .map(|i| ((w >> (i * LIMB_BITS)) & &self.limb_mask))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }

    pub fn find_w_modulus_ceil(&self, times: u64) -> [N; LIMBS] {
        let max = &self.w_ceil * (times + 1);
        let (n, rem) = max.div_rem(&self.w_modulus);
        let n = if rem.gt(&BigUint::from(0u64)) {
            n + 1u64
        } else {
            n
        };

        let mut upper = &self.w_modulus * n;

        let mut limbs = vec![];
        for _ in 0..LIMBS - 1 {
            let rem = (&upper & &self.limb_mask) + &self.limb_modulus * times;
            upper = (upper - &rem) >> LIMB_BITS;
            limbs.push(bn_to_field::<N>(&rem));
        }
        limbs.push(bn_to_field::<N>(&upper));
        limbs.try_into().unwrap()
    }
}
