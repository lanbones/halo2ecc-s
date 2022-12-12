use halo2_proofs::arithmetic::{BaseExt, FieldExt};
use num_bigint::BigUint;
use num_integer::Integer;

use super::{base_chip::BaseChipOps, range_chip::RangeChipOps};
use crate::{
    assign::{AssignedCondition, AssignedInteger, AssignedValue},
    context::Context,
    pair,
    utils::{bn_to_field, field_to_bn},
};

pub trait IntegerChipOps<W: BaseExt, N: FieldExt> {
    fn base_chip(&self) -> &dyn BaseChipOps<N>;
    fn range_chip(&self) -> &dyn RangeChipOps<W, N>;
    fn assign_w(&mut self, w: &BigUint) -> AssignedInteger<W, N>;
    fn assign_d(&mut self, v: &BigUint) -> (Vec<AssignedValue<N>>, AssignedValue<N>);
    fn conditionally_reduce(&mut self, a: AssignedInteger<W, N>) -> AssignedInteger<W, N>;
    fn reduce(&mut self, a: &AssignedInteger<W, N>) -> AssignedInteger<W, N>;
    fn int_add(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N>;
    fn int_sub(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N>;
    fn int_neg(&mut self, a: &AssignedInteger<W, N>) -> AssignedInteger<W, N>;
    fn int_mul(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N>;
    fn int_unsafe_invert(&mut self, x: &AssignedInteger<W, N>) -> AssignedInteger<W, N>;
    fn int_div(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> (AssignedCondition<N>, AssignedInteger<W, N>);
    fn is_pure_zero(&mut self, a: &AssignedInteger<W, N>) -> AssignedCondition<N>;
    fn is_pure_w_modulus(&mut self, a: &AssignedInteger<W, N>) -> AssignedCondition<N>;
    fn is_int_zero(&mut self, a: &AssignedInteger<W, N>) -> AssignedCondition<N>;
    fn is_int_equal(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedCondition<N> {
        let diff = self.int_sub(a, b);
        self.is_int_zero(&diff)
    }
    fn assign_int_constant(&mut self, w: W) -> AssignedInteger<W, N>;
    fn assert_int_equal(&mut self, a: &AssignedInteger<W, N>, b: &AssignedInteger<W, N>);
    fn int_square(&mut self, a: &AssignedInteger<W, N>) -> AssignedInteger<W, N>;
    fn int_mul_small_constant(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: u64,
    ) -> AssignedInteger<W, N>;
    fn bisec_int(
        &mut self,
        cond: &AssignedCondition<N>,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N>;
    fn get_w(&self, a: &AssignedInteger<W, N>) -> W;
}

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    fn add_constraints_for_mul_equation_on_limbs(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
        d: &Vec<AssignedValue<N>>,
        rem: &AssignedInteger<W, N>,
    ) {
        assert!(a.times < self.info().overflow_limit);
        assert!(b.times < self.info().overflow_limit);
        assert!(rem.times < self.info().overflow_limit);

        let info = self.info();
        let one = N::one();

        let mut limbs = vec![];
        for pos in 0..self.info().limbs as usize {
            // e.g. l0 = a0 * b0 - d0 * w0
            // e.g. l1 = a1 * b0 + a0 * b1 - d1 * w0 - d0 * w1
            // ...
            let l = self.mul_add_with_next_line(
                (0..pos + 1)
                    .map(|i| {
                        (
                            &a.limbs_le[i],
                            &b.limbs_le[pos - i],
                            &d[i],
                            info.neg_w_modulus_limbs_le[pos - i],
                        )
                    })
                    .collect(),
            );

            limbs.push(l);
        }

        // 4 self.info().limbs can reduce loop by merging 2 limbs check.

        // check sum limb[0]
        let u = self.sum_with_constant(
            vec![(&limbs[0], one), (&rem.limbs_le[0], -one)],
            Some(info.limb_modulus_n),
        );
        let (v, r) = field_to_bn(&u.val).div_rem(&info.limb_modulus);
        assert_eq!(r, BigUint::from(0u64));
        let (v_h_bn, v_l_bn) = v.div_rem(&info.limb_modulus);
        let mut v_h = self.assign_common(&v_h_bn);
        let mut v_l = self.assign_nonleading_limb(&v_l_bn);

        self.one_line_with_last(
            vec![
                pair!(&v_h, info.limb_coeffs[2]),
                pair!(&v_l, info.limb_coeffs[1]),
            ],
            pair!(&u, -one),
            None,
            (vec![], None),
        );

        // check sum limb[1..] with carry
        for i in 1..self.info().limbs as usize {
            let u = self.sum_with_constant(
                vec![(&limbs[i], one), (&rem.limbs_le[i], -one)],
                Some(info.limb_modulus_n - one),
            );
            let u = self.sum_with_constant(
                vec![
                    (&v_h, info.limb_coeffs[1]),
                    (&v_l, info.limb_coeffs[0]),
                    (&u, one),
                ],
                None,
            );

            let (v, r) = field_to_bn(&u.val).div_rem(&info.limb_modulus);
            assert_eq!(r, BigUint::from(0u64));
            let (v_h_bn, v_l_bn) = v.div_rem(&info.limb_modulus);
            v_h = self.assign_common(&v_h_bn);
            v_l = self.assign_nonleading_limb(&v_l_bn);

            BaseChipOps::one_line_with_last(
                self,
                vec![
                    pair!(&v_h, info.limb_coeffs[2]),
                    pair!(&v_l, info.limb_coeffs[1]),
                ],
                pair!(&u, -one),
                None,
                (vec![], None),
            );
        }
    }

    fn add_constraints_for_mul_equation_on_native(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
        d_native: &AssignedValue<N>,
        rem: &AssignedInteger<W, N>,
    ) {
        let info = self.info();
        let zero = N::zero();
        let one = N::one();
        self.one_line(
            vec![
                pair!(&a.native, zero),
                pair!(&b.native, zero),
                pair!(d_native, info.w_native),
                pair!(&rem.native, one),
            ],
            None,
            (vec![-one], None),
        );
    }

    fn get_w_bn(&self, a: &AssignedInteger<W, N>) -> BigUint {
        let mut res = BigUint::from(0u64);
        for i in (0..self.info().limbs as usize).rev() {
            res = res << self.info().limb_bits;
            res = res + field_to_bn(&a.limbs_le[i].val)
        }
        res
    }
}

impl<W: BaseExt, N: FieldExt> IntegerChipOps<W, N> for Context<W, N> {
    fn base_chip(&self) -> &dyn BaseChipOps<N> {
        self
    }

    fn range_chip(&self) -> &dyn RangeChipOps<W, N> {
        self
    }

    fn assign_w(&mut self, w: &BigUint) -> AssignedInteger<W, N> {
        let info = self.info();

        let mut limbs = vec![];
        for i in 0..self.info().limbs as u64 - 1 {
            limbs.push(
                self.assign_nonleading_limb(
                    &((w >> (i * self.info().limb_bits)) & &info.limb_mask),
                ),
            );
        }
        limbs.push(self.assign_w_ceil_leading_limb(
            &(w >> ((self.info().limbs as u64 - 1) * self.info().limb_bits) & &info.limb_mask),
        ));

        let schemas = limbs.iter().zip(info.limb_coeffs.clone());
        let native = self.sum_with_constant(schemas.collect(), None);

        AssignedInteger::new(limbs.try_into().unwrap(), native, 1)
    }

    fn assign_d(&mut self, d: &BigUint) -> (Vec<AssignedValue<N>>, AssignedValue<N>) {
        let info = self.info();

        let mut limbs = vec![];
        for i in 0..self.info().limbs as u64 - 1 {
            limbs.push(
                self.assign_nonleading_limb(
                    &((d >> (i * self.info().limb_bits)) & &info.limb_mask),
                ),
            );
        }
        limbs.push(self.assign_d_leading_limb(
            &(d >> ((self.info().limbs as u64 - 1) * self.info().limb_bits) & &info.limb_mask),
        ));

        let schemas = limbs.iter().zip(info.limb_coeffs.clone());
        let native = self.sum_with_constant(schemas.collect(), None);
        (limbs.try_into().unwrap(), native)
    }

    fn reduce(&mut self, a: &AssignedInteger<W, N>) -> AssignedInteger<W, N> {
        if a.times == 1 {
            return a.clone();
        }

        assert!(a.times < self.info().overflow_limit);

        let info = self.info();
        let one = N::one();

        // Check a = d * w + rem
        let a_bn = self.get_w_bn(&a);
        let (d, rem) = a_bn.div_rem(&info.w_modulus);

        // self.info().overflow_limit * limb_modulus > a.times * limb_modulus > a[0]
        // u = d * w[0] + rem[0] + self.info().overflow_limit * limb_modulus - a[0]
        //   < common * limb_modulus + limb_modulus + self.info().overflow_limit * limb_modulus
        // v = u / limbs < common
        //   < common + 1 + self.info().overflow_limit
        let u = &d * &info.w_modulus_limbs_le_bn[0]
            + &info.bn_to_limb_le(&rem)[0]
            + &info.limb_modulus * self.info().overflow_limit
            - field_to_bn(&a.limbs_le[0].val);

        let (v, v_rem) = u.div_rem(&info.limb_modulus);
        assert!(v_rem == BigUint::from(0u64));

        let rem = self.assign_w(&rem);
        let d = self.assign_common(&d);
        let v = self.assign_nonleading_limb(&v);

        // Constrain on native.
        self.one_line_with_last(
            vec![pair!(&d, info.w_native), pair!(&rem.native, one)],
            pair!(&a.native, -one),
            None,
            (vec![], None),
        );

        // constrains on limb_modulus
        self.one_line_with_last(
            vec![
                pair!(&d, info.w_modulus_limbs_le[0]),
                pair!(&rem.limbs_le[0], one),
                pair!(&a.limbs_le[0], -one),
            ],
            pair!(&v, -bn_to_field::<N>(&info.limb_modulus)),
            Some(bn_to_field(
                &(&info.limb_modulus * self.info().overflow_limit),
            )),
            (vec![], None),
        );

        rem
    }

    fn conditionally_reduce(&mut self, a: AssignedInteger<W, N>) -> AssignedInteger<W, N> {
        let threshold = 1 << (self.info().overflow_bits - 2);
        if a.times > threshold {
            self.reduce(&a)
        } else {
            a
        }
    }

    fn int_add(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N> {
        let info = self.info();
        let mut limbs = vec![];

        for i in 0..self.info().limbs as usize {
            let value = self.add(&a.limbs_le[i], &b.limbs_le[i]);
            limbs.push(value)
        }

        let schemas = limbs.iter().zip(info.limb_coeffs.clone());
        let native = self.sum_with_constant(schemas.collect(), None);

        let res = AssignedInteger::new(limbs.try_into().unwrap(), native, a.times + b.times);

        self.conditionally_reduce(res)
    }

    fn int_sub(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N> {
        let info = self.info();
        let upper_limbs = self.info().find_w_modulus_ceil(b.times);

        let one = N::one();

        let mut limbs = vec![];
        for i in 0..self.info().limbs as usize {
            let cell = self.sum_with_constant(
                vec![(&a.limbs_le[i], one), (&b.limbs_le[i], -one)],
                Some(upper_limbs[i]),
            );
            limbs.push(cell);
        }

        let schemas = limbs.iter().zip(info.limb_coeffs.clone());
        let native = self.sum_with_constant(schemas.collect(), None);

        let res = AssignedInteger::new(limbs.try_into().unwrap(), native, a.times + b.times + 2);
        self.conditionally_reduce(res)
    }

    fn int_neg(&mut self, a: &AssignedInteger<W, N>) -> AssignedInteger<W, N> {
        let info = self.info();
        let upper_limbs = self.info().find_w_modulus_ceil(a.times);

        let one = N::one();

        let mut limbs = vec![];
        for i in 0..self.info().limbs as usize {
            let cell = self.sum_with_constant(vec![(&a.limbs_le[i], -one)], Some(upper_limbs[i]));
            limbs.push(cell);
        }

        let schemas = limbs.iter().zip(info.limb_coeffs.clone());
        let native = self.sum_with_constant(schemas.collect(), None);

        let res = AssignedInteger::new(limbs.try_into().unwrap(), native, a.times + 2);
        self.conditionally_reduce(res)
    }

    fn int_mul(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N> {
        let info = self.info();
        let a_bn = self.get_w_bn(&a);
        let b_bn = self.get_w_bn(&b);
        let (d, rem) = (a_bn * b_bn).div_rem(&info.w_modulus);

        let rem = self.assign_w(&rem);
        let d = self.assign_d(&d);

        self.add_constraints_for_mul_equation_on_limbs(a, b, &d.0, &rem);
        self.add_constraints_for_mul_equation_on_native(a, b, &d.1, &rem);

        rem
    }
    fn int_unsafe_invert(&mut self, x: &AssignedInteger<W, N>) -> AssignedInteger<W, N> {
        //TODO: optimize
        let one = self.assign_int_constant(W::one());
        let (c, v) = self.int_div(&one, x);
        self.assert_false(&c);
        v
    }
    fn int_div(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> (AssignedCondition<N>, AssignedInteger<W, N>) {
        let info = self.info();

        // If b != 0
        // Find (c, d) that b * c = d * w + reduced_a,
        // Call reduce on `a` because if b = 1, we cannot find such (c, d), c < w_ceil and d >= 0

        // If b == 0
        // Find (c, d) that b * c = d * w + reduced_a * 0,

        let b = self.reduce(b);
        let is_b_zero = self.is_int_zero(&b);
        let a_coeff = self.not(&is_b_zero);

        let a = {
            let a = self.reduce(a);
            let mut limbs_le = vec![];
            for i in 0..self.info().limbs as usize {
                let cell = self.mul(&a.limbs_le[i], &a_coeff.0);
                limbs_le.push(cell);
            }
            let native = self.mul(&a.native, &a_coeff.0);
            AssignedInteger::<W, N>::new(limbs_le.try_into().unwrap(), native, a.times)
        };

        let a_bn = self.get_w_bn(&a);
        let b_bn = self.get_w_bn(&b);
        let c = bn_to_field::<W>(&b_bn)
            .invert()
            .map(|b| bn_to_field::<W>(&a_bn) * b)
            .unwrap_or(W::zero());
        let c_bn = field_to_bn(&c);
        let d_bn = (&b_bn * &c_bn - &a_bn) / &info.w_modulus;

        let c = self.assign_w(&c_bn);
        let d = self.assign_d(&d_bn);

        self.add_constraints_for_mul_equation_on_limbs(&b, &c, &d.0, &a);
        self.add_constraints_for_mul_equation_on_native(&b, &c, &d.1, &a);

        (is_b_zero, c)
    }

    fn is_pure_zero(&mut self, a: &AssignedInteger<W, N>) -> AssignedCondition<N> {
        let one = N::one();

        let sum = self.sum_with_constant(a.limbs_le.iter().map(|v| (v, one)).collect(), None);
        self.is_zero(&sum)
    }

    fn is_pure_w_modulus(&mut self, a: &AssignedInteger<W, N>) -> AssignedCondition<N> {
        let info = self.info();

        let native_diff = self.add_constant(&a.native, -info.w_native);
        let is_native_eq = self.is_zero(&native_diff);

        let limb0_diff = self.add_constant(&a.limbs_le[0], -info.w_modulus_limbs_le[0]);
        let is_limb0_eq = self.is_zero(&limb0_diff);

        self.and(&is_native_eq, &is_limb0_eq)
    }

    fn is_int_zero(&mut self, a: &AssignedInteger<W, N>) -> AssignedCondition<N> {
        let a = self.reduce(a);
        let is_zero = self.is_pure_zero(&a);
        let is_w_modulus = self.is_pure_w_modulus(&a);

        self.or(&is_zero, &is_w_modulus)
    }

    fn assign_int_constant(&mut self, w: W) -> AssignedInteger<W, N> {
        let info = self.info();

        let w = field_to_bn(&w);
        let limbs_value = info.bn_to_limb_le_n(&w);

        let mut limbs = vec![];
        for limb in limbs_value {
            let cell = self.assign_constant(limb);
            limbs.push(cell);
        }

        let native = self.assign_constant(bn_to_field(&(w % &info.n_modulus)));

        AssignedInteger::new(limbs.try_into().unwrap(), native, 1)
    }

    fn assert_int_equal(&mut self, a: &AssignedInteger<W, N>, b: &AssignedInteger<W, N>) {
        let zero = N::zero();

        let diff = self.int_sub(a, b);
        let diff = self.reduce(&diff);

        self.assert_constant(&diff.native, zero);
        self.assert_constant(&diff.limbs_le[0], zero);
    }

    fn int_square(&mut self, a: &AssignedInteger<W, N>) -> AssignedInteger<W, N> {
        self.int_mul(a, a)
    }

    fn int_mul_small_constant(
        &mut self,
        a: &AssignedInteger<W, N>,
        b: u64,
    ) -> AssignedInteger<W, N> {
        let threshold = 1 << (self.info().overflow_bits - 2);
        assert!(b < threshold);

        let info = self.info();

        let a_opt = if a.times * b >= self.info().overflow_limit {
            Some(self.reduce(a))
        } else {
            None
        };

        let a = if let Some(reduced_a) = &a_opt {
            reduced_a
        } else {
            a
        };

        let mut limbs = vec![];
        for i in 0..self.info().limbs as usize {
            let cell = self.sum_with_constant(vec![(&a.limbs_le[i], N::from(b as u64))], None);
            limbs.push(cell);
        }

        let schemas = limbs.iter().zip(info.limb_coeffs.clone());
        let native = self.sum_with_constant(schemas.collect(), None);

        let res = AssignedInteger::new(limbs.try_into().unwrap(), native, a.times * b);
        let res = self.conditionally_reduce(res);
        res
    }

    fn bisec_int(
        &mut self,
        cond: &AssignedCondition<N>,
        a: &AssignedInteger<W, N>,
        b: &AssignedInteger<W, N>,
    ) -> AssignedInteger<W, N> {
        let mut limbs = vec![];
        for i in 0..self.info().limbs as usize {
            let cell = self.bisec(cond, &a.limbs_le[i], &b.limbs_le[i]);
            limbs.push(cell);
        }
        let native = self.bisec(cond, &a.native, &b.native);

        AssignedInteger::new(
            limbs.try_into().unwrap(),
            native,
            u64::max(a.times, b.times),
        )
    }

    fn get_w(&self, a: &AssignedInteger<W, N>) -> W {
        bn_to_field(&self.get_w_bn(a))
    }
}
