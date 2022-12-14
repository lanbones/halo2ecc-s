use halo2_proofs::arithmetic::{BaseExt, FieldExt};
use num_bigint::BigUint;

use crate::{
    assign::{AssignedFq12, AssignedFq2, AssignedFq6},
    context::IntegerContext,
    utils::bn_to_field,
};

use super::{bn256_constants::{
    FROBENIUS_COEFF_FQ12_C1, FROBENIUS_COEFF_FQ2_C1, FROBENIUS_COEFF_FQ6_C1, FROBENIUS_COEFF_FQ6_C2,
}, integer_chip::IntegerChipOps};

impl<W: BaseExt, N: FieldExt> IntegerContext<W, N> {
    pub fn fq2_assert_equal(&mut self, x: &AssignedFq2<W, N>, y: &AssignedFq2<W, N>) {
        self.assert_int_equal(&x.0, &y.0);
        self.assert_int_equal(&x.1, &y.1);
    }
    pub fn fq2_assign_zero(&mut self) -> AssignedFq2<W, N> {
        let fq2_zero = self.assign_int_constant(W::zero());
        (fq2_zero.clone(), fq2_zero)
    }
    pub fn fq2_assign_one(&mut self) -> AssignedFq2<W, N> {
        (
            self.assign_int_constant(W::one()),
            self.assign_int_constant(W::zero()),
        )
    }
    pub fn fq2_assign_constant(&mut self, c0: W, c1: W) -> AssignedFq2<W, N> {
        (self.assign_int_constant(c0), self.assign_int_constant(c1))
    }
    pub fn fq2_add(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_add(&a.0, &b.0), self.int_add(&a.1, &b.1))
    }
    pub fn fq2_mul(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let ab00 = self.int_mul(&a.0, &b.0);
        let ab11 = self.int_mul(&a.1, &b.1);
        let c0 = self.int_sub(&ab00, &ab11);

        let a01 = self.int_add(&a.0, &a.1);
        let b01 = self.int_add(&b.0, &b.1);
        let c1 = self.int_mul(&a01, &b01);
        let c1 = self.int_sub(&c1, &ab00);
        let c1 = self.int_sub(&c1, &ab11);

        (c0, c1)
    }
    pub fn fq2_sub(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_sub(&a.0, &b.0), self.int_sub(&a.1, &b.1))
    }
    pub fn fq2_double(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_add(&a.0, &a.0), self.int_add(&a.1, &a.1))
    }
    pub fn fq2_square(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        self.fq2_mul(a, a)
    }
    pub fn fq2_neg(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_neg(&a.0), self.int_neg(&a.1))
    }
    pub fn fq2_conjugate(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (a.0.clone(), self.int_neg(&a.1))
    }
    pub fn fq2_mul_by_nonresidue(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let a2 = self.fq2_double(a);
        let a4 = self.fq2_double(&a2);
        let a8 = self.fq2_double(&a4);

        let t = self.int_add(&a8.0, &a.0);
        let c0 = self.int_sub(&t, &a.1);

        let t = self.int_add(&a8.1, &a.0);
        let c1 = self.int_add(&t, &a.1);

        (c0, c1)
    }
    pub fn fq2_unsafe_invert(&mut self, x: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let t0 = self.int_square(&x.0);
        let t1 = self.int_square(&x.1);
        let t0 = self.int_add(&t0, &t1);
        let t = self.int_unsafe_invert(&t0);
        let c0 = self.int_mul(&x.0, &t);
        let c1 = self.int_mul(&x.1, &t);
        let c1 = self.int_neg(&c1);
        (c0, c1)
    }
    pub fn fq2_frobenius_map(&mut self, x: &AssignedFq2<W, N>, power: usize) -> AssignedFq2<W, N> {
        let v = self.assign_int_constant(bn_to_field(&BigUint::from_bytes_le(
            &FROBENIUS_COEFF_FQ2_C1[power % 2],
        )));
        (x.0.clone(), self.int_mul(&x.1, &v))
    }
}

impl<W: BaseExt, N: FieldExt> IntegerContext<W, N> {
    pub fn fq6_assert_equal(&mut self, x: &AssignedFq6<W, N>, y: &AssignedFq6<W, N>) {
        self.fq2_assert_equal(&x.0, &y.0);
        self.fq2_assert_equal(&x.1, &y.1);
        self.fq2_assert_equal(&x.2, &y.2);
    }
    pub fn fq6_assign_zero(&mut self) -> AssignedFq6<W, N> {
        let fq2_zero = self.fq2_assign_zero();
        (fq2_zero.clone(), fq2_zero.clone(), fq2_zero)
    }
    pub fn fq6_assign_one(&mut self) -> AssignedFq6<W, N> {
        let fq2_one = self.fq2_assign_one();
        let fq2_zero = self.fq2_assign_zero();
        (fq2_one, fq2_zero.clone(), fq2_zero)
    }
    pub fn fq6_add(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_add(&a.0, &b.0),
            self.fq2_add(&a.1, &b.1),
            self.fq2_add(&a.2, &b.2),
        )
    }
    pub fn fq6_mul(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        let ab00 = self.fq2_mul(&a.0, &b.0);
        let ab11 = self.fq2_mul(&a.1, &b.1);
        let ab22 = self.fq2_mul(&a.2, &b.2);

        let c0 = {
            let b12 = self.fq2_add(&b.1, &b.2);
            let a12 = self.fq2_add(&a.1, &a.2);
            let t = self.fq2_mul(&a12, &b12);
            let t = self.fq2_sub(&t, &ab11);
            let t = self.fq2_sub(&t, &ab22);
            let t = self.fq2_mul_by_nonresidue(&t);
            self.fq2_add(&t, &ab00)
        };

        let c1 = {
            let b01 = self.fq2_add(&b.0, &b.1);
            let a01 = self.fq2_add(&a.0, &a.1);
            let t = self.fq2_mul(&a01, &b01);
            let t = self.fq2_sub(&t, &ab00);
            let t = self.fq2_sub(&t, &ab11);
            let ab22 = self.fq2_mul_by_nonresidue(&ab22);
            self.fq2_add(&t, &ab22)
        };

        let c2 = {
            let b02 = self.fq2_add(&b.0, &b.2);
            let a02 = self.fq2_add(&a.0, &a.2);
            let t = self.fq2_mul(&a02, &b02);
            let t = self.fq2_sub(&t, &ab00);
            let t = self.fq2_add(&t, &ab11);
            self.fq2_sub(&t, &ab22)
        };

        (c0, c1, c2)
    }
    pub fn fq6_sub(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_sub(&a.0, &b.0),
            self.fq2_sub(&a.1, &b.1),
            self.fq2_sub(&a.2, &b.2),
        )
    }
    pub fn fq6_double(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_double(&a.0),
            self.fq2_double(&a.1),
            self.fq2_double(&a.2),
        )
    }
    pub fn fq6_square(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        self.fq6_mul(a, a)
    }
    pub fn fq6_neg(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (self.fq2_neg(&a.0), self.fq2_neg(&a.1), self.fq2_neg(&a.2))
    }
    pub fn fq6_mul_by_nonresidue(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (self.fq2_mul_by_nonresidue(&a.2), a.0.clone(), a.1.clone())
    }
    pub fn fq6_mul_by_1(
        &mut self,
        a: &AssignedFq6<W, N>,
        b1: &AssignedFq2<W, N>,
    ) -> AssignedFq6<W, N> {
        let ab11 = self.fq2_mul(&a.1, &b1);

        let c0 = {
            let b12 = b1;
            let a12 = self.fq2_add(&a.1, &a.2);
            let t = self.fq2_mul(&a12, &b12);
            let t = self.fq2_sub(&t, &ab11);
            self.fq2_mul_by_nonresidue(&t)
        };

        let c1 = {
            let b01 = b1;
            let a01 = self.fq2_add(&a.0, &a.1);
            let t = self.fq2_mul(&a01, &b01);
            self.fq2_sub(&t, &ab11)
        };

        let c2 = ab11;

        (c0, c1, c2)
    }
    pub fn fq6_mul_by_01(
        &mut self,
        a: &AssignedFq6<W, N>,
        b0: &AssignedFq2<W, N>,
        b1: &AssignedFq2<W, N>,
    ) -> AssignedFq6<W, N> {
        let ab00 = self.fq2_mul(&a.0, &b0);
        let ab11 = self.fq2_mul(&a.1, &b1);

        let c0 = {
            let b12 = b1;
            let a12 = self.fq2_add(&a.1, &a.2);
            let t = self.fq2_mul(&a12, &b12);
            let t = self.fq2_sub(&t, &ab11);
            let t = self.fq2_mul_by_nonresidue(&t);
            self.fq2_add(&t, &ab00)
        };

        let c1 = {
            let b01 = self.fq2_add(b0, b1);
            let a01 = self.fq2_add(&a.0, &a.1);
            let t = self.fq2_mul(&a01, &b01);
            let t = self.fq2_sub(&t, &ab00);
            self.fq2_sub(&t, &ab11)
        };

        let c2 = {
            let b02 = b0;
            let a02 = self.fq2_add(&a.0, &a.2);
            let t = self.fq2_mul(&a02, &b02);
            let t = self.fq2_sub(&t, &ab00);
            self.fq2_add(&t, &ab11)
        };

        (c0, c1, c2)
    }
    pub fn fq6_unsafe_invert(&mut self, x: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        let c0 = self.fq2_mul_by_nonresidue(&x.2);
        let c0 = self.fq2_mul(&c0, &x.1);
        let c0 = self.fq2_neg(&c0);
        let x0s = self.fq2_square(&x.0);
        let c0 = self.fq2_add(&c0, &x0s);

        let c1 = self.fq2_square(&x.2);
        let c1 = self.fq2_mul_by_nonresidue(&c1);
        let x01 = self.fq2_mul(&x.0, &x.1);
        let c1 = self.fq2_sub(&c1, &x01);

        let c2 = self.fq2_square(&x.1);
        let x02 = self.fq2_mul(&x.0, &x.2);
        let c2 = self.fq2_sub(&c2, &x02);

        let c0x0 = self.fq2_mul(&c0, &x.0);
        let c1x2 = self.fq2_mul(&c1, &x.2);
        let c2x1 = self.fq2_mul(&c2, &x.1);
        let t = self.fq2_add(&c1x2, &c2x1);
        let t = self.fq2_mul_by_nonresidue(&t);
        let t = self.fq2_add(&t, &c0x0);

        let t = self.fq2_unsafe_invert(&t);

        (
            self.fq2_mul(&t, &c0),
            self.fq2_mul(&t, &c1),
            self.fq2_mul(&t, &c2),
        )
    }
    pub fn fq6_frobenius_map(&mut self, x: &AssignedFq6<W, N>, power: usize) -> AssignedFq6<W, N> {
        let c0 = self.fq2_frobenius_map(&x.0, power);
        let c1 = self.fq2_frobenius_map(&x.1, power);
        let c2 = self.fq2_frobenius_map(&x.2, power);

        let coeff_c1 =
            FROBENIUS_COEFF_FQ6_C1[power % 6].map(|x| bn_to_field(&BigUint::from_bytes_le(&x)));
        let coeff_c1 = self.fq2_assign_constant(coeff_c1[0], coeff_c1[1]);
        let c1 = self.fq2_mul(&c1, &coeff_c1);
        let coeff_c2 =
            FROBENIUS_COEFF_FQ6_C2[power % 6].map(|x| bn_to_field(&BigUint::from_bytes_le(&x)));
        let coeff_c2 = self.fq2_assign_constant(coeff_c2[0], coeff_c2[1]);
        let c2 = self.fq2_mul(&c2, &coeff_c2);

        (c0, c1, c2)
    }
}

impl<W: BaseExt, N: FieldExt> IntegerContext<W, N> {
    pub fn fq12_assert_one(&mut self, x: &AssignedFq12<W, N>) {
        let one = self.fq12_assign_one();
        self.fq12_assert_eq(x, &one);
    }
    pub fn fq12_assert_eq(&mut self, x: &AssignedFq12<W, N>, y: &AssignedFq12<W, N>) {
        self.fq6_assert_equal(&x.0, &y.0);
        self.fq6_assert_equal(&x.1, &y.1);
    }
    pub fn fq12_assign_zero(&mut self) -> AssignedFq12<W, N> {
        let fq6_zero = self.fq6_assign_zero();
        (fq6_zero.clone(), fq6_zero)
    }
    pub fn fq12_assign_one(&mut self) -> AssignedFq12<W, N> {
        let fq6_one = self.fq6_assign_one();
        let fq6_zero = self.fq6_assign_zero();
        (fq6_one, fq6_zero)
    }
    pub fn fq12_add(
        &mut self,
        a: &AssignedFq12<W, N>,
        b: &AssignedFq12<W, N>,
    ) -> AssignedFq12<W, N> {
        (self.fq6_add(&a.0, &b.0), self.fq6_add(&a.1, &b.1))
    }
    pub fn fq12_mul(
        &mut self,
        a: &AssignedFq12<W, N>,
        b: &AssignedFq12<W, N>,
    ) -> AssignedFq12<W, N> {
        let ab00 = self.fq6_mul(&a.0, &b.0);
        let ab11 = self.fq6_mul(&a.1, &b.1);

        let a01 = self.fq6_add(&a.0, &a.1);
        let b01 = self.fq6_add(&b.0, &b.1);
        let c1 = self.fq6_mul(&a01, &b01);
        let c1 = self.fq6_sub(&c1, &ab00);
        let c1 = self.fq6_sub(&c1, &ab11);

        let ab11 = self.fq6_mul_by_nonresidue(&ab11);
        let c0 = self.fq6_add(&ab00, &ab11);

        (c0, c1)
    }
    pub fn fq12_sub(
        &mut self,
        a: &AssignedFq12<W, N>,
        b: &AssignedFq12<W, N>,
    ) -> AssignedFq12<W, N> {
        (self.fq6_sub(&a.0, &b.0), self.fq6_sub(&a.1, &b.1))
    }
    pub fn fq12_double(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_double(&a.0), self.fq6_double(&a.1))
    }
    pub fn fq12_square(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        self.fq12_mul(a, a)
    }
    pub fn fq12_neg(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_neg(&a.0), self.fq6_neg(&a.1))
    }
    pub fn fq12_conjugate(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (x.0.clone(), self.fq6_neg(&x.1))
    }
    pub fn fq12_mul_by_014(
        &mut self,
        x: &AssignedFq12<W, N>,
        c0: &AssignedFq2<W, N>,
        c1: &AssignedFq2<W, N>,
        c4: &AssignedFq2<W, N>,
    ) -> AssignedFq12<W, N> {
        let t0 = self.fq6_mul_by_01(&x.0, c0, c1);
        let t1 = self.fq6_mul_by_1(&x.1, c4);
        let o = self.fq2_add(c1, c4);

        let x0 = self.fq6_add(&x.0, &x.1);
        let x0 = self.fq6_mul_by_01(&x0, c0, &o);
        let x0 = self.fq6_sub(&x0, &t0);
        let x0 = self.fq6_sub(&x0, &t1);

        let x1 = self.fq6_mul_by_nonresidue(&t1);
        let x1 = self.fq6_add(&x1, &t0);

        (x0, x1)
    }
    pub fn fq12_mul_by_034(
        &mut self,
        x: &AssignedFq12<W, N>,
        c0: &AssignedFq2<W, N>,
        c3: &AssignedFq2<W, N>,
        c4: &AssignedFq2<W, N>,
    ) -> AssignedFq12<W, N> {
        let t00 = self.fq2_mul(&x.0 .0, c0);
        let t01 = self.fq2_mul(&x.0 .1, c0);
        let t02 = self.fq2_mul(&x.0 .2, c0);
        let t0 = (t00, t01, t02);

        let t1 = self.fq6_mul_by_01(&x.1, c3, c4);
        let t2 = self.fq6_add(&x.0, &x.1);
        let o = self.fq2_add(c0, c3);
        let t2 = self.fq6_mul_by_01(&t2, &o, c4);
        let t2 = self.fq6_sub(&t2, &t0);
        let x1 = self.fq6_sub(&t2, &t1);
        let t1 = self.fq6_mul_by_nonresidue(&t1);
        let x0 = self.fq6_add(&t0, &t1);
        (x0, x1)
    }
    pub fn fp4_square(
        &mut self,
        c0: &mut AssignedFq2<W, N>,
        c1: &mut AssignedFq2<W, N>,
        a0: &AssignedFq2<W, N>,
        a1: &AssignedFq2<W, N>,
    ) {
        let t0 = self.fq2_square(&a0);
        let t1 = self.fq2_square(&a1);
        let mut t2 = self.fq2_mul_by_nonresidue(&t1);
        *c0 = self.fq2_add(&t2, &t0);
        t2 = self.fq2_add(a0, a1);
        t2 = self.fq2_square(&t2);
        t2 = self.fq2_sub(&t2, &t0);
        *c1 = self.fq2_sub(&t2, &t1);
    }
    pub fn fq12_cyclotomic_square(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        let zero = self.fq2_assign_zero();
        let mut t3 = zero.clone();
        let mut t4 = zero.clone();
        let mut t5 = zero.clone();
        let mut t6 = zero;

        self.fp4_square(&mut t3, &mut t4, &x.0 .0, &x.1 .1);
        let mut t2 = self.fq2_sub(&t3, &x.0 .0);
        t2 = self.fq2_double(&t2);
        let c00 = self.fq2_add(&t2, &t3);

        t2 = self.fq2_add(&t4, &x.1 .1);
        t2 = self.fq2_double(&t2);
        let c11 = self.fq2_add(&t2, &t4);

        self.fp4_square(&mut t3, &mut t4, &x.1 .0, &x.0 .2);
        self.fp4_square(&mut t5, &mut t6, &x.0 .1, &x.1 .2);

        t2 = self.fq2_sub(&t3, &x.0 .1);
        t2 = self.fq2_double(&t2);
        let c01 = self.fq2_add(&t2, &t3);
        t2 = self.fq2_add(&t4, &x.1 .2);
        t2 = self.fq2_double(&t2);
        let c12 = self.fq2_add(&t2, &t4);
        t3 = t6;
        t3 = self.fq2_mul_by_nonresidue(&t3);
        t2 = self.fq2_add(&t3, &x.1 .0);
        t2 = self.fq2_double(&t2);
        let c10 = self.fq2_add(&t2, &t3);
        t2 = self.fq2_sub(&t5, &x.0 .2);
        t2 = self.fq2_double(&t2);
        let c02 = self.fq2_add(&t2, &t5);

        ((c00, c01, c02), (c10, c11, c12))
    }
    pub fn fq12_unsafe_invert(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        let x0s = self.fq6_square(&x.0);
        let x1s = self.fq6_square(&x.1);
        let t = self.fq6_mul_by_nonresidue(&x1s);
        let t = self.fq6_sub(&x0s, &t);
        let t = self.fq6_unsafe_invert(&t);

        let c0 = self.fq6_mul(&t, &x.0);
        let c1 = self.fq6_mul(&t, &x.1);
        let c1 = self.fq6_neg(&c1);
        (c0, c1)
    }
    pub fn fq12_frobenius_map(
        &mut self,
        x: &AssignedFq12<W, N>,
        power: usize,
    ) -> AssignedFq12<W, N> {
        let c0 = self.fq6_frobenius_map(&x.0, power);
        let c1 = self.fq6_frobenius_map(&x.1, power);

        let coeff =
            FROBENIUS_COEFF_FQ12_C1[power % 12].map(|x| bn_to_field(&BigUint::from_bytes_le(&x)));
        let coeff = self.fq2_assign_constant(coeff[0], coeff[1]);
        let c1c0 = self.fq2_mul(&c1.0, &coeff);
        let c1c1 = self.fq2_mul(&c1.1, &coeff);
        let c1c2 = self.fq2_mul(&c1.2, &coeff);

        (c0, (c1c0, c1c1, c1c2))
    }
}
