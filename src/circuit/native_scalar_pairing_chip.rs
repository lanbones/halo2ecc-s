use std::marker::PhantomData;

use crate::assign::{AssignedFq6, AssignedG2, AssignedG2Affine, AssignedG2Prepared};
use crate::circuit::integer_chip::IntegerChipOps;
use crate::context::Context;
use crate::utils::bn_to_field;
use crate::{
    assign::{AssignedFq12, AssignedFq2, AssignedG1Affine},
    context::NativeScalarEccContext,
};
use halo2_proofs::arithmetic::{BaseExt, CurveAffine, FieldExt};
use num_bigint::BigUint;

use super::base_chip::BaseChipOps;

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    fn fq2_assign_one(&mut self) -> AssignedFq2<W, N> {
        (
            self.assign_int_constant(W::one()),
            self.assign_int_constant(W::zero()),
        )
    }
    fn fq2_assign_constant(&mut self, c0: W, c1: W) -> AssignedFq2<W, N> {
        (self.assign_int_constant(c0), self.assign_int_constant(c1))
    }
    fn fq2_add(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_add(&a.0, &b.0), self.int_add(&a.1, &b.1))
    }
    fn fq2_mul(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
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
    fn fq2_sub(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_sub(&a.0, &b.0), self.int_sub(&a.1, &b.1))
    }
    fn fq2_double(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_add(&a.0, &a.0), self.int_add(&a.1, &a.1))
    }
    fn fq2_square(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        self.fq2_mul(a, a)
    }
    fn fq2_neg(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_neg(&a.0), self.int_neg(&a.1))
    }
    fn fq2_conjugate(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (a.0, self.int_neg(&a.1))
    }
    fn fq2_mul_by_nonresidue(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let a2 = self.fq2_double(a);
        let a4 = self.fq2_double(&a2);
        let a8 = self.fq2_double(&a4);

        let t = self.int_add(&a8.0, &a.0);
        let c0 = self.int_sub(&t, &a.1);

        let t = self.int_add(&a8.1, &a.0);
        let c1 = self.int_add(&t, &a.1);

        (c0, c1)
    }
}

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    fn fq6_add(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_add(&a.0, &b.0),
            self.fq2_add(&a.1, &b.1),
            self.fq2_add(&a.2, &b.2),
        )
    }
    fn fq6_mul(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
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
            let t = self.fq2_mul_by_nonresidue(&t);
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
    fn fq6_sub(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_sub(&a.0, &b.0),
            self.fq2_sub(&a.1, &b.1),
            self.fq2_sub(&a.2, &b.2),
        )
    }
    fn fq6_double(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_double(&a.0),
            self.fq2_double(&a.1),
            self.fq2_double(&a.2),
        )
    }
    fn fq6_square(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        self.fq6_mul(a, a)
    }
    fn fq6_neg(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (self.fq2_neg(&a.0), self.fq2_neg(&a.1), self.fq2_neg(&a.2))
    }
    fn fq6_mul_by_nonresidue(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (self.fq2_mul_by_nonresidue(&a.2), a.0, a.1)
    }
    fn fq6_mul_by_1(&mut self, a: &AssignedFq6<W, N>, b1: &AssignedFq2<W, N>) -> AssignedFq6<W, N> {
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
    fn fq6_mul_by_01(
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
}

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    fn fq12_add(&mut self, a: &AssignedFq12<W, N>, b: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_add(&a.0, &b.0), self.fq6_add(&a.1, &b.1))
    }
    fn fq12_mul(&mut self, a: &AssignedFq12<W, N>, b: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
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
    fn fq12_sub(&mut self, a: &AssignedFq12<W, N>, b: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_sub(&a.0, &b.0), self.fq6_sub(&a.1, &b.1))
    }
    fn fq12_double(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_double(&a.0), self.fq6_double(&a.1))
    }
    fn fq12_square(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        self.fq12_mul(a, a)
    }
    fn fq12_neg(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_neg(&a.0), self.fq6_neg(&a.1))
    }
    fn fq12_conjugate(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (x.0, self.fq6_neg(&x.1))
    }
    fn fq12_mul_by_014(
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
    fn fq12_mul_by_034(
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
        let x0 = self.fq6_sub(&t2, &t1);
        let t1 = self.fq6_mul_by_nonresidue(&t1);
        let x1 = self.fq6_add(&t0, &t1);
        (x0, x1)
    }
}

impl<C: CurveAffine> NativeScalarEccContext<C> {
    fn doubling_step(
        &mut self,
        pt: &AssignedG2<C, C::Scalar>,
    ) -> (
        AssignedG2<C, C::Scalar>,
        [AssignedFq2<C::Base, C::Scalar>; 3],
    ) {
        //see https://eprint.iacr.org/2010/354.pdf
        let x2 = self.0.fq2_square(&pt.x);

        let y2 = self.0.fq2_square(&pt.y);
        let _2y2 = self.0.fq2_double(&y2);
        let _4y2 = self.0.fq2_double(&_2y2);
        let _4y4 = self.0.fq2_square(&_2y2);
        let _8y4 = self.0.fq2_double(&_4y4);

        let z2 = self.0.fq2_square(&pt.z);

        let _4xy2 = {
            let t = self.0.fq2_mul(&y2, &pt.x);
            let t = self.0.fq2_double(&t);
            let t = self.0.fq2_double(&t);
            t
        };
        let _3x2 = {
            let t = self.0.fq2_double(&x2);
            let t = self.0.fq2_add(&t, &x2);
            t
        };
        let _6x2 = self.0.fq2_double(&_3x2);
        let _9x4 = self.0.fq2_square(&_3x2);
        let _3x2_x = self.0.fq2_add(&_3x2, &pt.x);

        let rx = {
            let t = self.0.fq2_sub(&_9x4, &_4xy2);
            let t = self.0.fq2_sub(&t, &_4xy2);
            t
        };

        let ry = {
            let t = self.0.fq2_sub(&_4xy2, &pt.x);
            let t = self.0.fq2_mul(&t, &_3x2_x);
            let t = self.0.fq2_sub(&t, &_8y4);
            t
        };

        let rz = {
            let yz = self.0.fq2_mul(&pt.y, &pt.z);
            self.0.fq2_double(&yz)
        };

        let c0 = {
            let t = self.0.fq2_mul(&z2, &rz);
            self.0.fq2_double(&t)
        };

        let c1 = {
            let _6x2z2 = self.0.fq2_mul(&z2, &_6x2);
            self.0.fq2_neg(&_6x2z2)
        };

        let c2 = {
            let _6x3 = self.0.fq2_mul(&_6x2, &pt.x);
            self.0.fq2_sub(&_6x3, &_4y2)
        };

        (AssignedG2::new(rx, ry, rz), [c0, c1, c2])
    }

    fn addition_step(
        &mut self,
        pt: &AssignedG2<C, C::Scalar>,
        pq: &AssignedG2Affine<C, C::Scalar>,
    ) -> (
        AssignedG2<C, C::Scalar>,
        [AssignedFq2<C::Base, C::Scalar>; 3],
    ) {
        //see https://eprint.iacr.org/2010/354.pdf
        let zt2 = self.0.fq2_square(&pt.z);
        let yqzt = self.0.fq2_mul(&pq.y, &pt.z);
        let yqzt3 = self.0.fq2_mul(&yqzt, &zt2);
        let _2yqzt3 = self.0.fq2_double(&yqzt3);
        let _2yt = self.0.fq2_double(&pt.y);
        let _2yqzt3_2yt = self.0.fq2_sub(&_2yqzt3, &_2yt);

        let xqzt2 = self.0.fq2_mul(&pq.x, &zt2);
        let xqzt2_xt = self.0.fq2_sub(&xqzt2, &pt.x);
        let _2_xqzt2_xt = self.0.fq2_double(&xqzt2_xt); // 2(xqzt2 - xt)
        let _4_xqzt2_xt_2 = self.0.fq2_square(&_2_xqzt2_xt); // 4(xqzt2 - xt) ^ 2

        let rx = {
            let t0 = self.0.fq2_mul(&_4_xqzt2_xt_2, &xqzt2_xt); // 4(xqzt2 - xt) ^ 3
            let t1 = self.0.fq2_double(&_4_xqzt2_xt_2); // 8(xqzt2 - xt) ^ 2
            let t2 = self.0.fq2_mul(&t1, &pt.x); // 8(xqzt2 - xt) ^ 2 * xt

            let t = self.0.fq2_square(&_2yqzt3_2yt); // (2yqzt3 - 2yt) ^ 2
            let t = self.0.fq2_sub(&t, &t0); // (2yqzt3 - 2yt) ^ 2 - 4(xqzt2 - xt) ^ 3
            let t = self.0.fq2_sub(&t, &t2); // (2yqzt3 - 2yt) ^ 2 - 4(xqzt2 - xt) ^ 3 - 8(xqzt2 - xt) ^ 2 * xt
            t
        };

        let ry = {
            let t0 = self.0.fq2_mul(&_4_xqzt2_xt_2, &pt.x); // 4(xqzt2 - xt) ^ 2 * xt
            let t0 = self.0.fq2_sub(&t0, &rx); // 4(xqzt2 - xt) ^ 2 * xt - xr
            let t0 = self.0.fq2_mul(&_2yqzt3_2yt, &t0); // (2yqzt3 - 2yt) * (4(xqzt2 - xt) ^ 2 * xt - xr)
            let t1 = self.0.fq2_mul(&_2_xqzt2_xt, &_4_xqzt2_xt_2); // 8(xqzt2 - xt) ^ 3
            let t1 = self.0.fq2_mul(&t1, &pt.y); // 8yt * (xqzt2 - xt) ^ 3
            let t = self.0.fq2_sub(&t0, &t1);
            t
        };

        let rz = self.0.fq2_mul(&pt.z, &_2_xqzt2_xt);

        let c0 = self.0.fq2_double(&rz);
        let c1 = {
            let t = self.0.fq2_add(&_2yqzt3, &_2yt); // 2(yqzt3 + yt)
            self.0.fq2_double(&t)
        };
        let c2 = {
            let t0 = self.0.fq2_mul(&_2yqzt3, &pq.x); // 2yqzt3xq
            let t0 = self.0.fq2_sub(&t0, &_2yt); // 2(yqzt3xq - yt)
            let t0 = self.0.fq2_double(&t0); // 4(yqzt3xq - yt)
            let t0 = self.0.fq2_mul(&t0, &pq.x); // 4xq(yqzt3xq - yt)
            let t1 = self.0.fq2_mul(&pq.y, &rz); // yqzr
            let t1 = self.0.fq2_double(&t1); // 2yqzr
            self.0.fq2_sub(&t0, &t1) // 4xq(yqzt3xq - yt) - yqzr
        };

        (AssignedG2::new(rx, ry, rz), [c0, c1, c2])
    }

    fn g2affine_to_g2(&mut self, g2: &AssignedG2Affine<C, C::Scalar>) -> AssignedG2<C, C::Scalar> {
        // not support identity
        self.0.assert_false(&g2.z);
        let z = self.0.fq2_assign_one();

        AssignedG2::new(g2.x, g2.y, z)
    }

    fn g2_neg(&mut self, g2: &AssignedG2Affine<C, C::Scalar>) -> AssignedG2Affine<C, C::Scalar> {
        let y = self.0.fq2_neg(&g2.y);
        AssignedG2Affine::new(g2.x, y, g2.z)
    }

    fn prepare_g2(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2Prepared<C, C::ScalarExt> {
        // currently only supports bn256
        // 6U+2 for in NAF form
        pub const SIX_U_PLUS_2_NAF: [i8; 65] = [
            0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0,
            0, 0, 1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1,
            1, 0, -1, 0, 0, 1, 0, 1, 1,
        ];

        let neg_g2 = self.g2_neg(&g2);

        let mut coeffs = vec![];
        let mut r = self.g2affine_to_g2(g2);

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            let (new_r, coeff) = self.doubling_step(&r);
            r = new_r;
            coeffs.push(coeff);
            let x = SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    let (new_r, coeff) = self.addition_step(&r, &g2);
                    r = new_r;
                    coeffs.push(coeff);
                }
                -1 => {
                    let (new_r, coeff) = self.addition_step(&r, &neg_g2);
                    r = new_r;
                    coeffs.push(coeff);
                }
                _ => continue,
            }
        }

        let mut q1 = g2.clone();

        // TODO: fill constant value
        let FROBENIUS_COEFF_FQ6_C1_1 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(&[])),
            bn_to_field(&BigUint::from_bytes_le(&[])),
        );let FROBENIUS_COEFF_FQ6_C1_2 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(&[])),
            bn_to_field(&BigUint::from_bytes_le(&[])),
        );
        let XI_TO_Q_MINUS_1_OVER_2 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(&[])),
            bn_to_field(&BigUint::from_bytes_le(&[])),
        );

        q1.x.1 = self.0.int_neg(&q1.x.1);
        q1.x = self.0.fq2_mul(&q1.x, &FROBENIUS_COEFF_FQ6_C1_1);

        q1.y.1 = self.0.int_neg(&q1.y.1);
        q1.y = self.0.fq2_mul(&q1.y, &XI_TO_Q_MINUS_1_OVER_2);

        let (new_r, coeff) = self.addition_step(&r, &q1);
        r = new_r;
        coeffs.push(coeff);

        let mut minusq2 = g2.clone();
        minusq2.x = self.0.fq2_mul(&minusq2.x, &FROBENIUS_COEFF_FQ6_C1_2);

        let (new_r, coeff) = self.addition_step(&r, &minusq2);
        r = new_r;
        coeffs.push(coeff);

        AssignedG2Prepared::new(coeffs)
    }

    fn ell(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
        coeffs: &(
            AssignedFq2<C::Base, C::Scalar>,
            AssignedFq2<C::Base, C::Scalar>,
            AssignedFq2<C::Base, C::Scalar>,
        ),
        p: &AssignedG1Affine<C, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let c00 = &coeffs.0 .0;
        let c01 = &coeffs.0 .1;
        let c10 = &coeffs.1 .0;
        let c11 = &coeffs.1 .1;

        let c00 = self.0.int_mul(&c00, &p.y);
        let c01 = self.0.int_mul(&c01, &p.y);
        let c10 = self.0.int_mul(&c10, &p.x);
        let c11 = self.0.int_mul(&c11, &p.x);

        self.0
            .fq12_mul_by_034(f, &(c00, c01), &(c10, c11), &coeffs.2)
    }
}
