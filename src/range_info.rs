use halo2_proofs::arithmetic::BaseExt;
use halo2_proofs::arithmetic::FieldExt;
use num_bigint::BigUint;
use num_integer::Integer;
use std::marker::PhantomData;

use crate::circuit::range_chip::COMMON_RANGE_BITS;
use crate::circuit::range_chip::MAX_CHUNKS;
use crate::utils::{bn_to_field, field_to_bn};

#[derive(Debug, Clone)]
pub struct RangeInfo<W: BaseExt, N: FieldExt> {
    pub limbs: u64,
    pub limb_bits: u64,

    pub w_ceil_leading_chunks: u64,
    pub n_floor_leading_chunks: u64,
    pub d_leading_chunks: u64,

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

    pub w_modulus_limbs_le_bn: Vec<BigUint>,
    pub w_modulus_limbs_le: Vec<N>,
    pub neg_w_modulus_limbs_le: Vec<N>,
    pub limb_coeffs: Vec<N>,
    pub limb_modulus_n: N,

    pub overflow_bits: u64,
    pub overflow_limit: u64,

    pub w_native: N,

    pub pure_w_check_limbs: u64,
    pub reduce_check_limbs: u64,
    pub w_modulus_of_ceil_times: Vec<Option<Vec<N>>>,

    pub _phantom: PhantomData<W>,
}

impl<W: BaseExt, N: FieldExt> RangeInfo<W, N> {
    fn bits_to_leading_bits_and_chunks(bits: u64, common_bits: u64) -> (u64, u64) {
        let leading_chunk_bits = bits % common_bits;
        if leading_chunk_bits == 0 {
            (common_bits, bits / common_bits % MAX_CHUNKS)
        } else {
            (leading_chunk_bits, bits / common_bits % MAX_CHUNKS + 1)
        }
    }

    pub fn new(common_bits: u64, overflow_bits: u64) -> Self {
        let w_max = field_to_bn(&-W::one());
        let w_ceil_bits = w_max.bits();
        assert!(BigUint::from(1u64) << w_ceil_bits > w_max);
        assert!(BigUint::from(1u64) << (w_ceil_bits - 1) < w_max);
        let (w_ceil_leading_bits, w_ceil_leading_chunks) =
            Self::bits_to_leading_bits_and_chunks(w_ceil_bits, common_bits);

        let n_max = field_to_bn(&-N::one());
        let n_floor_bits = n_max.bits() - 1;
        assert!(BigUint::from(1u64) << n_floor_bits < n_max);
        assert!(BigUint::from(1u64) << (n_floor_bits + 1) >= n_max);
        let (n_floor_leading_bits, n_floor_leading_chunks) =
            Self::bits_to_leading_bits_and_chunks(n_floor_bits, common_bits);

        let d_bits = Self::d_bits(overflow_bits);
        let (d_leading_bits, d_leading_chunks) =
            Self::bits_to_leading_bits_and_chunks(d_bits, common_bits);

        let limb_bits = common_bits * MAX_CHUNKS;
        let limbs = (w_ceil_bits + limb_bits - 1) / limb_bits;

        let max_d = BigUint::from(1u64) << d_bits;
        let limb_mask = (BigUint::from(1u64) << limb_bits) - 1u64;
        let n_modulus = &n_max + 1u64;
        let w_modulus = &w_max + 1u64;
        let w_native = &w_modulus % &n_modulus;
        let integer_modulus = BigUint::from(1u64) << (limbs * limb_bits);
        let neg_w_modulus = &integer_modulus - &w_modulus;

        let w_modulus_limbs_le_bn = (0..limbs)
            .into_iter()
            .map(|i| ((&w_modulus >> (i * limb_bits)) & &limb_mask))
            .collect::<Vec<_>>();
        let neg_w_modulus_limbs_le = (0..limbs)
            .into_iter()
            .map(|i| (bn_to_field::<N>(&((&neg_w_modulus >> (i * limb_bits)) & &limb_mask))))
            .collect::<Vec<_>>();
        let w_modulus_limbs_le = w_modulus_limbs_le_bn
            .iter()
            .map(|x| bn_to_field(x))
            .collect::<Vec<_>>();

        let limb_modulus = BigUint::from(1u64) << limb_bits;
        let limb_modulus_n = bn_to_field(&limb_modulus);

        let overflow_limit = 1u64 << overflow_bits;
        let w_modulus_of_ceil_times = vec![None; overflow_limit as usize];

        let mut res = Self {
            w_ceil_bits,
            d_leading_bits,
            w_ceil_leading_bits,
            n_floor_leading_bits,

            max_d,
            w_ceil: BigUint::from(1u64) << w_ceil_bits,
            n_modulus,
            w_modulus,
            neg_w_modulus,
            common_range_mask: BigUint::from((1u64 << common_bits) - 1),
            limb_mask,
            limb_modulus,
            integer_modulus,
            limb_coeffs: (0..limbs)
                .into_iter()
                .map(|i| bn_to_field(&(BigUint::from(1u64) << (i * limb_bits))))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),

            limb_modulus_n,
            neg_w_modulus_limbs_le,
            w_modulus_limbs_le,
            w_modulus_limbs_le_bn,

            w_native: bn_to_field::<N>(&w_native),

            limbs,
            limb_bits,
            w_ceil_leading_chunks,
            n_floor_leading_chunks,
            d_leading_chunks,

            overflow_bits,
            overflow_limit,

            pure_w_check_limbs: (w_ceil_bits - n_floor_bits + limb_bits - 1) / limb_bits,
            reduce_check_limbs: (w_ceil_bits + overflow_bits - n_floor_bits + limb_bits - 1)
                / limb_bits,
            w_modulus_of_ceil_times,

            _phantom: PhantomData,
        };

        for i in 1..overflow_limit {
            res.w_modulus_of_ceil_times[i as usize] = Some(res.find_w_modulus_of_ceil_times(i));
        }

        res.pre_check();
        res
    }

    fn pre_check(&self) {
        // is_pure_w_modulus():
        // lcm(limb, native) >= w_ceil
        let limb_check_modulus = BigUint::from(1u64) << (self.limb_bits * self.pure_w_check_limbs);
        let lcm = self.n_modulus.lcm(&limb_check_modulus);
        assert!(lcm >= self.w_ceil);

        // reduce():
        // lcm(limb ^ reduce_check_limbs, native) >= w_ceil * self.overflow_limit
        let limb_modulus = BigUint::from(1u64) << (self.limb_bits * self.reduce_check_limbs);
        let lcm = self.n_modulus.lcm(&limb_modulus);
        assert!(lcm >= &self.w_ceil * self.overflow_limit);
        // Ensure that d can be assigned by assign_common.
        assert!(COMMON_RANGE_BITS > self.overflow_bits);
        // Ensure that v can be assigned by assign_nonleading_limb
        assert!(
            &((BigUint::from(1u64) << COMMON_RANGE_BITS) + 2u64 + self.overflow_bits)
                <= &self.limb_modulus
        );
        // Ensure that v * limb_modulus would not overflow
        assert!(&(&self.limb_modulus * &self.limb_modulus) < &self.n_modulus);

        // mul():
        // lcm(integer_modulus, native) >= w_ceil * w_ceil * self.overflow_limit * self.overflow_limit
        let lcm = self.integer_modulus.lcm(&self.n_modulus);
        let max_a = &self.w_modulus * self.overflow_limit;
        let max_b = &self.w_modulus * self.overflow_limit;
        let max_l = max_a * max_b;
        let max_d = &self.max_d;
        let max_w = &self.w_modulus;
        let max_rem = &self.w_ceil;
        let max_r = max_d * max_w + max_rem;
        assert!(max_l <= lcm);
        assert!(max_r <= lcm);
        assert!(max_l <= max_r);
        // (u, v)
        let sum_limb_max = self.limbs
            * (self.overflow_limit * self.overflow_limit + 1)
            * &self.limb_modulus
            * &self.limb_modulus;
        assert!(sum_limb_max < self.n_modulus);
        // Ensure we can split v into (common, nonleading)
        // The last `1` reserve for the borrow and carry.
        // The real reserved value is (1 + common) / limb_modulus.
        assert!(
            1 << COMMON_RANGE_BITS
                > self.limbs * (self.overflow_limit * self.overflow_limit + 1 + 1)
        );
        let common_modulus = BigUint::from(1u64) << COMMON_RANGE_BITS;
        assert!(&common_modulus + 1u64 < self.limb_modulus);

        // Ensure nonoverflow for one limb sum check.
        assert!(&self.limb_modulus * &self.limb_modulus * common_modulus < self.n_modulus);

        // algorithm limitation
        assert!(self.limbs >= 3);
    }

    pub fn d_bits(overflow_bits: u64) -> u64 {
        // a * b = w * d + rem

        // -- a <= w_ceil_bits + overflow_bits /\ b <= w_ceil_bits + overflow_bits
        // -> a * b <= w_ceil_bits * 2 + overflow_bits * 2 - 1

        // -- w * d + rem >= a * b
        // <- d_bits + w_ceil_bits - 2 = w_ceil_bits * 2 + overflow_bits * 2 - 1

        let w_max = field_to_bn(&-W::one());
        let w_ceil_bits = w_max.bits();
        let d_bits = w_ceil_bits + overflow_bits * 2 + 1;
        assert!(d_bits + w_ceil_bits - 2 >= w_ceil_bits * 2 + overflow_bits * 2 - 1);
        d_bits
    }

    pub fn bn_to_limb_le_n(&self, w: &BigUint) -> Vec<N> {
        (0..self.limbs)
            .into_iter()
            .map(|i| bn_to_field(&((w >> (i * self.limb_bits)) & &self.limb_mask)))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }

    pub fn bn_to_limb_le(&self, w: &BigUint) -> Vec<BigUint> {
        (0..self.limbs)
            .into_iter()
            .map(|i| ((w >> (i * self.limb_bits)) & &self.limb_mask))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }

    pub fn find_w_modulus_of_ceil_times(&self, times: u64) -> Vec<N> {
        let max = &self.w_ceil * times;
        let (n, rem) = max.div_rem(&self.w_modulus);
        let n = if rem.gt(&BigUint::from(0u64)) {
            n + 1u64
        } else {
            n
        };

        let mut upper = &self.w_modulus * n;

        let mut limbs = vec![];
        for _ in 0..self.limbs - 1 {
            let rem = (&upper & &self.limb_mask) + &self.limb_modulus * times;
            upper = (upper - &rem) >> self.limb_bits;
            limbs.push(bn_to_field::<N>(&rem));
        }

        assert!(upper >= (BigUint::from(1u64) << (self.w_ceil_bits % self.limb_bits)) * times);
        assert!(upper < (BigUint::from(1u64) << (self.w_ceil_bits % self.limb_bits)) * (times + 1));

        limbs.push(bn_to_field::<N>(&upper));
        limbs.try_into().unwrap()
    }
}

#[test]
fn test_range_info() {
    {
        use halo2_proofs::pairing::bn256::Fq;
        use halo2_proofs::pairing::bn256::Fr;

        RangeInfo::<Fq, Fr>::new(18, 6);
    }

    {
        use halo2_proofs::pairing::bls12_381::Fr as Bls12_381_Fr;
        use halo2_proofs::pairing::bn256::Fr;

        RangeInfo::<Bls12_381_Fr, Fr>::new(18, 6);
        //println!("info {:?}", info);
    }

    {
        //use halo2_proofs::pairing::bls12_381::Fq as Bls12_381_Fq;
        //use halo2_proofs::pairing::bn256::Fr;

        //let info = RangeInfo::<Bls12_381_Fq, Fr>::new(18, 6);
        //println!("info {:?}", info);
    }
}
