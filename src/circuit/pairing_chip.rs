/*
  The implementation is ported from https://github.com/privacy-scaling-explorations/pairing
*/
use crate::assign::{AssignedFq12, AssignedFq2, AssignedG1Affine};
use crate::assign::{AssignedG2, AssignedG2Affine, AssignedG2Prepared};
use halo2_proofs::arithmetic::{CurveAffine, FieldExt};

use super::fq12::{Fq12BnSpecificOps, Fq12ChipOps};

pub trait PairingChipOps<C: CurveAffine, N: FieldExt>:
    Fq12ChipOps<C::Base, N> + Fq12BnSpecificOps<C::Base, N>
{
    fn doubling_step(&mut self, pt: &mut AssignedG2<C, N>) -> [AssignedFq2<C::Base, N>; 3] {
        //see https://eprint.iacr.org/2010/354.pdf
        let x2 = self.fq2_square(&pt.x);

        let y2 = self.fq2_square(&pt.y);
        let _2y2 = self.fq2_double(&y2);
        let _4y2 = self.fq2_double(&_2y2);
        let _4y4 = self.fq2_square(&_2y2);
        let _8y4 = self.fq2_double(&_4y4);

        let z2 = self.fq2_square(&pt.z);

        let _4xy2 = {
            let t = self.fq2_mul(&y2, &pt.x);
            let t = self.fq2_double(&t);
            let t = self.fq2_double(&t);
            t
        };
        let _3x2 = {
            let t = self.fq2_double(&x2);
            let t = self.fq2_add(&t, &x2);
            t
        };
        let _6x2 = self.fq2_double(&_3x2);
        let _9x4 = self.fq2_square(&_3x2);
        let _3x2_x = self.fq2_add(&_3x2, &pt.x);

        let rx = {
            let t = self.fq2_sub(&_9x4, &_4xy2);
            let t = self.fq2_sub(&t, &_4xy2);
            t
        };

        let ry = {
            let t = self.fq2_sub(&_4xy2, &rx);
            let t = self.fq2_mul(&t, &_3x2);
            let t = self.fq2_sub(&t, &_8y4);
            t
        };

        let rz = {
            let yz = self.fq2_mul(&pt.y, &pt.z);
            self.fq2_double(&yz)
        };

        let c0 = {
            let t = self.fq2_mul(&z2, &rz);
            self.fq2_double(&t)
        };

        let c1 = {
            let _6x2z2 = self.fq2_mul(&z2, &_6x2);
            self.fq2_neg(&_6x2z2)
        };

        let c2 = {
            let _6x3 = self.fq2_mul(&_6x2, &pt.x);
            self.fq2_sub(&_6x3, &_4y2)
        };

        *pt = AssignedG2::new(rx, ry, rz);

        [c0, c1, c2]
    }

    fn addition_step(
        &mut self,
        pt: &mut AssignedG2<C, N>,
        pq: &AssignedG2Affine<C, N>,
    ) -> [AssignedFq2<C::Base, N>; 3] {
        //see https://eprint.iacr.org/2010/354.pdf
        let zt2 = self.fq2_square(&pt.z);
        let yqzt = self.fq2_mul(&pq.y, &pt.z);
        let yqzt3 = self.fq2_mul(&yqzt, &zt2);
        let yqzt3_yt = self.fq2_sub(&yqzt3, &pt.y);
        let _2yqzt3_2yt = self.fq2_double(&yqzt3_yt);

        let xqzt2 = self.fq2_mul(&pq.x, &zt2);
        let xqzt2_xt = self.fq2_sub(&xqzt2, &pt.x);
        let _2_xqzt2_xt = self.fq2_double(&xqzt2_xt); // 2(xqzt2 - xt)
        let _4_xqzt2_xt_2 = self.fq2_square(&_2_xqzt2_xt); // 4(xqzt2 - xt) ^ 2

        let rx = {
            let t0 = self.fq2_mul(&_4_xqzt2_xt_2, &xqzt2_xt); // 4(xqzt2 - xt) ^ 3
            let t1 = self.fq2_double(&_4_xqzt2_xt_2); // 8(xqzt2 - xt) ^ 2
            let t2 = self.fq2_mul(&t1, &pt.x); // 8(xqzt2 - xt) ^ 2 * xt

            let t = self.fq2_square(&_2yqzt3_2yt); // (2yqzt3 - 2yt) ^ 2
            let t = self.fq2_sub(&t, &t0); // (2yqzt3 - 2yt) ^ 2 - 4(xqzt2 - xt) ^ 3
            let t = self.fq2_sub(&t, &t2); // (2yqzt3 - 2yt) ^ 2 - 4(xqzt2 - xt) ^ 3 - 8(xqzt2 - xt) ^ 2 * xt
            t
        };

        let ry = {
            let t0 = self.fq2_mul(&_4_xqzt2_xt_2, &pt.x); // 4(xqzt2 - xt) ^ 2 * xt
            let t0 = self.fq2_sub(&t0, &rx); // 4(xqzt2 - xt) ^ 2 * xt - xr
            let t0 = self.fq2_mul(&_2yqzt3_2yt, &t0); // (2yqzt3 - 2yt) * (4(xqzt2 - xt) ^ 2 * xt - xr)
            let t1 = self.fq2_mul(&_2_xqzt2_xt, &_4_xqzt2_xt_2); // 8(xqzt2 - xt) ^ 3
            let t1 = self.fq2_mul(&t1, &pt.y); // 8yt * (xqzt2 - xt) ^ 3
            let t = self.fq2_sub(&t0, &t1);
            t
        };

        let rz = self.fq2_mul(&pt.z, &_2_xqzt2_xt);

        let c0 = self.fq2_double(&rz);
        let c1 = {
            let t = self.fq2_double(&_2yqzt3_2yt);
            self.fq2_neg(&t)
        };
        let c2 = {
            let t0 = self.fq2_double(&_2yqzt3_2yt); // 4(yqzt3xq - yt)
            let t0 = self.fq2_mul(&t0, &pq.x); // 4xq(yqzt3xq - yt)
            let t1 = self.fq2_mul(&pq.y, &rz); // yqzr
            let t1 = self.fq2_double(&t1); // 2yqzr
            self.fq2_sub(&t0, &t1) // 4xq(yqzt3xq - yt) - yqzr
        };

        *pt = AssignedG2::new(rx, ry, rz);
        [c0, c1, c2]
    }

    fn g2affine_to_g2(&mut self, g2: &AssignedG2Affine<C, N>) -> AssignedG2<C, N> {
        // not support identity
        self.base_integer_chip().base_chip().assert_false(&g2.z);
        let z = self.fq2_assign_one();

        AssignedG2::new(g2.x.clone(), g2.y.clone(), z)
    }

    fn g2_neg(&mut self, g2: &AssignedG2Affine<C, N>) -> AssignedG2Affine<C, N> {
        let y = self.fq2_neg(&g2.y);
        AssignedG2Affine::new(g2.x.clone(), y, g2.z)
    }

    fn prepare_g2(&mut self, g2: &AssignedG2Affine<C, N>) -> AssignedG2Prepared<C, N>;

    fn multi_miller_loop(
        &mut self,
        terms: &[(&AssignedG1Affine<C, N>, &AssignedG2Prepared<C, N>)],
    ) -> AssignedFq12<C::Base, N>;

    fn final_exponentiation(&mut self, f: &AssignedFq12<C::Base, N>) -> AssignedFq12<C::Base, N>;

    fn pairing(
        &mut self,
        terms: &[(&AssignedG1Affine<C, N>, &AssignedG2Affine<C, N>)],
    ) -> AssignedFq12<C::Base, N> {
        let prepared_terms = terms
            .iter()
            .map(|(p, q)| (*p, self.prepare_g2(q)))
            .collect::<Vec<_>>();
        let terms = prepared_terms
            .iter()
            .map(|(p, q)| (*p, q))
            .collect::<Vec<_>>();
        let res = self.multi_miller_loop(&terms[..]);
        self.final_exponentiation(&res)
    }

    fn check_pairing(&mut self, terms: &[(&AssignedG1Affine<C, N>, &AssignedG2Affine<C, N>)]) {
        let res = self.pairing(terms);
        self.fq12_assert_one(&res);
    }
}
