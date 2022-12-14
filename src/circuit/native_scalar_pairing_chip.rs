/*
  The implementation is ported from https://github.com/privacy-scaling-explorations/pairing
*/
use crate::assign::{AssignedG2, AssignedG2Affine, AssignedG2Prepared};
use crate::circuit::integer_chip::IntegerChipOps;
use crate::utils::bn_to_field;
use crate::{
    assign::{AssignedFq12, AssignedFq2, AssignedG1Affine},
    context::NativeScalarEccContext,
};
use halo2_proofs::arithmetic::CurveAffine;
use num_bigint::BigUint;

use super::base_chip::BaseChipOps;
use super::bn256_constants::*;

impl<C: CurveAffine> NativeScalarEccContext<C> {
    pub fn doubling_step(
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
            let t = self.0.fq2_sub(&_4xy2, &rx);
            let t = self.0.fq2_mul(&t, &_3x2);
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

    pub fn addition_step(
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
        let yqzt3_yt = self.0.fq2_sub(&yqzt3, &pt.y);
        let _2yqzt3_2yt = self.0.fq2_double(&yqzt3_yt);

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
            let t = self.0.fq2_double(&_2yqzt3_2yt);
            self.0.fq2_neg(&t)
        };
        let c2 = {
            let t0 = self.0.fq2_double(&_2yqzt3_2yt); // 4(yqzt3xq - yt)
            let t0 = self.0.fq2_mul(&t0, &pq.x); // 4xq(yqzt3xq - yt)
            let t1 = self.0.fq2_mul(&pq.y, &rz); // yqzr
            let t1 = self.0.fq2_double(&t1); // 2yqzr
            self.0.fq2_sub(&t0, &t1) // 4xq(yqzt3xq - yt) - yqzr
        };

        (AssignedG2::new(rx, ry, rz), [c0, c1, c2])
    }

    pub fn g2affine_to_g2(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2<C, C::Scalar> {
        // not support identity
        self.0.ctx.borrow_mut().assert_false(&g2.z);
        let z = self.0.fq2_assign_one();

        AssignedG2::new(g2.x.clone(), g2.y.clone(), z)
    }

    pub fn g2_neg(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2Affine<C, C::Scalar> {
        let y = self.0.fq2_neg(&g2.y);
        AssignedG2Affine::new(g2.x.clone(), y, g2.z)
    }

    pub fn prepare_g2(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2Prepared<C, C::ScalarExt> {
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

        let c11 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[1][0])),
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[1][1])),
        );
        let c12 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[2][0])),
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[2][1])),
        );
        let xi = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(&XI_TO_Q_MINUS_1_OVER_2[0])),
            bn_to_field(&BigUint::from_bytes_le(&XI_TO_Q_MINUS_1_OVER_2[1])),
        );

        q1.x.1 = self.0.int_neg(&q1.x.1);
        q1.x = self.0.fq2_mul(&q1.x, &c11);

        q1.y.1 = self.0.int_neg(&q1.y.1);
        q1.y = self.0.fq2_mul(&q1.y, &xi);

        let (new_r, coeff) = self.addition_step(&r, &q1);
        r = new_r;
        coeffs.push(coeff);

        let mut minusq2 = g2.clone();
        minusq2.x = self.0.fq2_mul(&minusq2.x, &c12);

        let (_, coeff) = self.addition_step(&r, &minusq2);
        coeffs.push(coeff);

        AssignedG2Prepared::new(coeffs)
    }

    pub fn ell(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
        coeffs: &[AssignedFq2<C::Base, C::Scalar>; 3],
        p: &AssignedG1Affine<C, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let c00 = &coeffs[0].0;
        let c01 = &coeffs[0].1;
        let c10 = &coeffs[1].0;
        let c11 = &coeffs[1].1;

        let c00 = self.0.int_mul(&c00, &p.y);
        let c01 = self.0.int_mul(&c01, &p.y);
        let c10 = self.0.int_mul(&c10, &p.x);
        let c11 = self.0.int_mul(&c11, &p.x);

        self.0
            .fq12_mul_by_034(f, &(c00, c01), &(c10, c11), &coeffs[2])
    }

    pub fn multi_miller_loop(
        &mut self,
        terms: &[(
            &AssignedG1Affine<C, C::Scalar>,
            &AssignedG2Prepared<C, C::Scalar>,
        )],
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let mut pairs = vec![];
        for &(p, q) in terms {
            // not support identity
            self.0.ctx.borrow_mut().assert_false(&p.z);
            pairs.push((p, q.coeffs.iter()));
        }

        let mut f = self.0.fq12_assign_one();

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            if i != SIX_U_PLUS_2_NAF.len() - 1 {
                f = self.0.fq12_square(&f);
            }
            for &mut (p, ref mut coeffs) in &mut pairs {
                f = self.ell(&f, coeffs.next().unwrap(), &p);
            }
            let x = SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    for &mut (p, ref mut coeffs) in &mut pairs {
                        f = self.ell(&f, coeffs.next().unwrap(), &p);
                    }
                }
                -1 => {
                    for &mut (p, ref mut coeffs) in &mut pairs {
                        f = self.ell(&f, coeffs.next().unwrap(), &p);
                    }
                }
                _ => continue,
            }
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            f = self.ell(&f, coeffs.next().unwrap(), &p);
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            f = self.ell(&f, coeffs.next().unwrap(), &p);
        }

        for &mut (_p, ref mut coeffs) in &mut pairs {
            assert!(coeffs.next().is_none());
        }

        f
    }

    pub fn exp_by_x(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let x = BN_X;
        let mut res = self.0.fq12_assign_one();
        for i in (0..64).rev() {
            res = self.0.fq12_cyclotomic_square(&res);
            if ((x >> i) & 1) == 1 {
                res = self.0.fq12_mul(&res, &f);
            }
        }
        res
    }

    pub fn final_exponentiation(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let f1 = self.0.fq12_conjugate(&f);
        let mut f2 = self.0.fq12_unsafe_invert(&f);

        let mut r = self.0.fq12_mul(&f1, &f2);
        f2 = r.clone();
        r = self.0.fq12_frobenius_map(&r, 2);
        r = self.0.fq12_mul(&r, &f2);

        let mut fp = r.clone();
        fp = self.0.fq12_frobenius_map(&fp, 1);

        let mut fp2 = r.clone();
        fp2 = self.0.fq12_frobenius_map(&fp2, 2);
        let mut fp3 = fp2.clone();
        fp3 = self.0.fq12_frobenius_map(&fp3, 1);

        let mut fu = r.clone();
        fu = self.exp_by_x(&fu);

        let mut fu2 = fu.clone();
        fu2 = self.exp_by_x(&fu2);

        let mut fu3 = fu2.clone();
        fu3 = self.exp_by_x(&fu3);

        let mut y3 = fu.clone();
        y3 = self.0.fq12_frobenius_map(&y3, 1);

        let mut fu2p = fu2.clone();
        fu2p = self.0.fq12_frobenius_map(&fu2p, 1);

        let mut fu3p = fu3.clone();
        fu3p = self.0.fq12_frobenius_map(&fu3p, 1);

        let mut y2 = fu2.clone();
        y2 = self.0.fq12_frobenius_map(&y2, 2);

        let mut y0 = fp;
        y0 = self.0.fq12_mul(&y0, &fp2);
        y0 = self.0.fq12_mul(&y0, &fp3);

        let mut y1 = r;
        y1 = self.0.fq12_conjugate(&y1);

        let mut y5 = fu2;
        y5 = self.0.fq12_conjugate(&y5);

        y3 = self.0.fq12_conjugate(&y3);

        let mut y4 = fu;
        y4 = self.0.fq12_mul(&y4, &fu2p);
        y4 = self.0.fq12_conjugate(&y4);

        let mut y6 = fu3;
        y6 = self.0.fq12_mul(&y6, &fu3p);
        y6 = self.0.fq12_conjugate(&y6);

        y6 = self.0.fq12_cyclotomic_square(&y6);
        y6 = self.0.fq12_mul(&y6, &y4);
        y6 = self.0.fq12_mul(&y6, &y5);

        let mut t1 = y3;
        t1 = self.0.fq12_mul(&t1, &y5);
        t1 = self.0.fq12_mul(&t1, &y6);

        y6 = self.0.fq12_mul(&y6, &y2);

        t1 = self.0.fq12_cyclotomic_square(&t1);
        t1 = self.0.fq12_mul(&t1, &y6);
        t1 = self.0.fq12_cyclotomic_square(&t1);

        let mut t0 = t1.clone();
        t0 = self.0.fq12_mul(&t0, &y1);

        t1 = self.0.fq12_mul(&t1, &y0);

        t0 = self.0.fq12_cyclotomic_square(&t0);
        t0 = self.0.fq12_mul(&t0, &t1);

        t0
    }

    pub fn check_pairing(
        &mut self,
        terms: &[(
            &AssignedG1Affine<C, C::Scalar>,
            &AssignedG2Affine<C, C::Scalar>,
        )],
    ) {
        let prepared_terms = terms
            .iter()
            .map(|(p, q)| (*p, self.prepare_g2(q)))
            .collect::<Vec<_>>();
        let terms = prepared_terms
            .iter()
            .map(|(p, q)| (*p, q))
            .collect::<Vec<_>>();
        let res = self.multi_miller_loop(&terms[..]);
        let res = self.final_exponentiation(&res);
        self.0.fq12_assert_one(&res);
    }
}
