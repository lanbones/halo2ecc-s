/*
  The implementation is ported from https://github.com/privacy-scaling-explorations/pairing
*/

use halo2_proofs::{
    arithmetic::FieldExt,
    pairing::bn256::{Fq, Fr, G1Affine},
};
use num_bigint::BigUint;

use crate::{
    assign::{
        AssignedFq12, AssignedFq2, AssignedFq6, AssignedG1Affine, AssignedG2Affine,
        AssignedG2Prepared,
    },
    context::NativeScalarEccContext,
    utils::bn_to_field,
};

use super::{
    bn256_constants::*,
    ecc_chip::EccBaseIntegerChipWrapper,
    fq12::{
        Fq12BnSpecificOps, Fq12ChipOps, Fq2BnSpecificOps, Fq2ChipOps, Fq6BnSpecificOps, Fq6ChipOps,
    },
    pairing_chip::PairingChipOps,
};

impl<N: FieldExt, T: EccBaseIntegerChipWrapper<Fq, N> + Fq2ChipOps<Fq, N>> Fq2BnSpecificOps<Fq, N>
    for T
{
    fn fq2_mul_by_nonresidue(&mut self, a: &AssignedFq2<Fq, N>) -> AssignedFq2<Fq, N> {
        let a2 = self.fq2_double(a);
        let a4 = self.fq2_double(&a2);
        let a8 = self.fq2_double(&a4);

        let t = self.base_integer_chip().int_add(&a8.0, &a.0);
        let c0 = self.base_integer_chip().int_sub(&t, &a.1);

        let t = self.base_integer_chip().int_add(&a8.1, &a.0);
        let c1 = self.base_integer_chip().int_add(&t, &a.1);

        (c0, c1)
    }

    fn fq2_frobenius_map(&mut self, x: &AssignedFq2<Fq, N>, power: usize) -> AssignedFq2<Fq, N> {
        let v = self
            .base_integer_chip()
            .assign_int_constant(bn_to_field(&BigUint::from_bytes_le(
                &FROBENIUS_COEFF_FQ2_C1[power % 2],
            )));
        (x.0.clone(), self.base_integer_chip().int_mul(&x.1, &v))
    }
}

impl<N: FieldExt, T: EccBaseIntegerChipWrapper<Fq, N> + Fq2ChipOps<Fq, N>> Fq6BnSpecificOps<Fq, N>
    for T
{
    fn fq6_mul_by_nonresidue(&mut self, a: &AssignedFq6<Fq, N>) -> AssignedFq6<Fq, N> {
        (self.fq2_mul_by_nonresidue(&a.2), a.0.clone(), a.1.clone())
    }

    fn fq6_frobenius_map(&mut self, x: &AssignedFq6<Fq, N>, power: usize) -> AssignedFq6<Fq, N> {
        let c0 = self.fq2_frobenius_map(&x.0, power);
        let c1 = self.fq2_frobenius_map(&x.1, power);
        let c2 = self.fq2_frobenius_map(&x.2, power);

        let coeff_c1 =
            FROBENIUS_COEFF_FQ6_C1[power % 6].map(|x| bn_to_field(&BigUint::from_bytes_le(&x)));
        let coeff_c1 = self.fq2_assign_constant((coeff_c1[0], coeff_c1[1]));
        let c1 = self.fq2_mul(&c1, &coeff_c1);
        let coeff_c2 =
            FROBENIUS_COEFF_FQ6_C2[power % 6].map(|x| bn_to_field(&BigUint::from_bytes_le(&x)));
        let coeff_c2 = self.fq2_assign_constant((coeff_c2[0], coeff_c2[1]));
        let c2 = self.fq2_mul(&c2, &coeff_c2);

        (c0, c1, c2)
    }
}

impl<N: FieldExt, T: EccBaseIntegerChipWrapper<Fq, N> + Fq2ChipOps<Fq, N> + Fq6ChipOps<Fq, N>>
    Fq12BnSpecificOps<Fq, N> for T
{
    fn fq12_frobenius_map(&mut self, x: &AssignedFq12<Fq, N>, power: usize) -> AssignedFq12<Fq, N> {
        let c0 = self.fq6_frobenius_map(&x.0, power);
        let c1 = self.fq6_frobenius_map(&x.1, power);

        let coeff =
            FROBENIUS_COEFF_FQ12_C1[power % 12].map(|x| bn_to_field(&BigUint::from_bytes_le(&x)));
        let coeff = self.fq2_assign_constant((coeff[0], coeff[1]));
        let c1c0 = self.fq2_mul(&c1.0, &coeff);
        let c1c1 = self.fq2_mul(&c1.1, &coeff);
        let c1c2 = self.fq2_mul(&c1.2, &coeff);

        (c0, (c1c0, c1c1, c1c2))
    }
}

impl Fq2ChipOps<Fq, Fr> for NativeScalarEccContext<G1Affine> {}
impl Fq6ChipOps<Fq, Fr> for NativeScalarEccContext<G1Affine> {}
impl Fq12ChipOps<Fq, Fr> for NativeScalarEccContext<G1Affine> {}

impl NativeScalarEccContext<G1Affine> {
    fn prepare_g2(
        &mut self,
        g2: &AssignedG2Affine<G1Affine, Fr>,
    ) -> AssignedG2Prepared<G1Affine, Fr> {
        let neg_g2 = self.g2_neg(&g2);

        let mut coeffs = vec![];
        let mut r = self.g2affine_to_g2(g2);

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            coeffs.push(self.doubling_step(&mut r));
            let x = SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    coeffs.push(self.addition_step(&mut r, &g2));
                }
                -1 => {
                    coeffs.push(self.addition_step(&mut r, &neg_g2));
                }
                _ => continue,
            }
        }

        let mut q1 = g2.clone();

        let c11 = self.fq2_assign_constant((
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[1][0])),
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[1][1])),
        ));
        let c12 = self.fq2_assign_constant((
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[2][0])),
            bn_to_field(&BigUint::from_bytes_le(&FROBENIUS_COEFF_FQ6_C1[2][1])),
        ));
        let xi = self.fq2_assign_constant((
            bn_to_field(&BigUint::from_bytes_le(&XI_TO_Q_MINUS_1_OVER_2[0])),
            bn_to_field(&BigUint::from_bytes_le(&XI_TO_Q_MINUS_1_OVER_2[1])),
        ));

        q1.x.1 = self.base_integer_chip().int_neg(&q1.x.1);
        q1.x = self.fq2_mul(&q1.x, &c11);

        q1.y.1 = self.base_integer_chip().int_neg(&q1.y.1);
        q1.y = self.fq2_mul(&q1.y, &xi);

        coeffs.push(self.addition_step(&mut r, &q1));

        let mut minusq2 = g2.clone();
        minusq2.x = self.fq2_mul(&minusq2.x, &c12);
        coeffs.push(self.addition_step(&mut r, &minusq2));

        AssignedG2Prepared::new(coeffs)
    }

    fn ell(
        &mut self,
        f: &AssignedFq12<Fq, Fr>,
        coeffs: &[AssignedFq2<Fq, Fr>; 3],
        p: &AssignedG1Affine<G1Affine, Fr>,
    ) -> AssignedFq12<Fq, Fr> {
        let c00 = &coeffs[0].0;
        let c01 = &coeffs[0].1;
        let c10 = &coeffs[1].0;
        let c11 = &coeffs[1].1;

        let c00 = self.base_integer_chip().int_mul(&c00, &p.y);
        let c01 = self.base_integer_chip().int_mul(&c01, &p.y);
        let c10 = self.base_integer_chip().int_mul(&c10, &p.x);
        let c11 = self.base_integer_chip().int_mul(&c11, &p.x);

        self.fq12_mul_by_034(f, &(c00, c01), &(c10, c11), &coeffs[2])
    }

    fn multi_miller_loop(
        &mut self,
        terms: &[(
            &AssignedG1Affine<G1Affine, Fr>,
            &AssignedG2Prepared<G1Affine, Fr>,
        )],
    ) -> AssignedFq12<Fq, Fr> {
        let mut pairs = vec![];
        for &(p, q) in terms {
            // not support identity
            self.base_integer_chip().base_chip().assert_false(&p.z);
            pairs.push((p, q.coeffs.iter()));
        }

        let mut f = self.fq12_assign_one();

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            if i != SIX_U_PLUS_2_NAF.len() - 1 {
                f = self.fq12_square(&f);
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

    fn exp_by_x(&mut self, f: &AssignedFq12<Fq, Fr>) -> AssignedFq12<Fq, Fr> {
        let x = BN_X;
        let mut res = self.fq12_assign_one();
        for i in (0..64).rev() {
            res = self.fq12_cyclotomic_square(&res);
            if ((x >> i) & 1) == 1 {
                res = self.fq12_mul(&res, &f);
            }
        }
        res
    }

    fn final_exponentiation(&mut self, f: &AssignedFq12<Fq, Fr>) -> AssignedFq12<Fq, Fr> {
        let f1 = self.fq12_conjugate(&f);
        let mut f2 = self.fq12_unsafe_invert(&f);

        let mut r = self.fq12_mul(&f1, &f2);
        f2 = r.clone();
        r = self.fq12_frobenius_map(&r, 2);
        r = self.fq12_mul(&r, &f2);

        let mut fp = r.clone();
        fp = self.fq12_frobenius_map(&fp, 1);

        let mut fp2 = r.clone();
        fp2 = self.fq12_frobenius_map(&fp2, 2);
        let mut fp3 = fp2.clone();
        fp3 = self.fq12_frobenius_map(&fp3, 1);

        let mut fu = r.clone();
        fu = self.exp_by_x(&fu);

        let mut fu2 = fu.clone();
        fu2 = self.exp_by_x(&fu2);

        let mut fu3 = fu2.clone();
        fu3 = self.exp_by_x(&fu3);

        let mut y3 = fu.clone();
        y3 = self.fq12_frobenius_map(&y3, 1);

        let mut fu2p = fu2.clone();
        fu2p = self.fq12_frobenius_map(&fu2p, 1);

        let mut fu3p = fu3.clone();
        fu3p = self.fq12_frobenius_map(&fu3p, 1);

        let mut y2 = fu2.clone();
        y2 = self.fq12_frobenius_map(&y2, 2);

        let mut y0 = fp;
        y0 = self.fq12_mul(&y0, &fp2);
        y0 = self.fq12_mul(&y0, &fp3);

        let mut y1 = r;
        y1 = self.fq12_conjugate(&y1);

        let mut y5 = fu2;
        y5 = self.fq12_conjugate(&y5);

        y3 = self.fq12_conjugate(&y3);

        let mut y4 = fu;
        y4 = self.fq12_mul(&y4, &fu2p);
        y4 = self.fq12_conjugate(&y4);

        let mut y6 = fu3;
        y6 = self.fq12_mul(&y6, &fu3p);
        y6 = self.fq12_conjugate(&y6);

        y6 = self.fq12_cyclotomic_square(&y6);
        y6 = self.fq12_mul(&y6, &y4);
        y6 = self.fq12_mul(&y6, &y5);

        let mut t1 = y3;
        t1 = self.fq12_mul(&t1, &y5);
        t1 = self.fq12_mul(&t1, &y6);

        y6 = self.fq12_mul(&y6, &y2);

        t1 = self.fq12_cyclotomic_square(&t1);
        t1 = self.fq12_mul(&t1, &y6);
        t1 = self.fq12_cyclotomic_square(&t1);

        let mut t0 = t1.clone();
        t0 = self.fq12_mul(&t0, &y1);

        t1 = self.fq12_mul(&t1, &y0);

        t0 = self.fq12_cyclotomic_square(&t0);
        t0 = self.fq12_mul(&t0, &t1);

        t0
    }
}

impl PairingChipOps<G1Affine, Fr> for NativeScalarEccContext<G1Affine> {
    fn prepare_g2(
        &mut self,
        g2: &AssignedG2Affine<G1Affine, Fr>,
    ) -> AssignedG2Prepared<G1Affine, Fr> {
        self.prepare_g2(g2)
    }

    fn multi_miller_loop(
        &mut self,
        terms: &[(
            &AssignedG1Affine<G1Affine, Fr>,
            &AssignedG2Prepared<G1Affine, Fr>,
        )],
    ) -> AssignedFq12<<G1Affine as halo2_proofs::arithmetic::CurveAffine>::Base, Fr> {
        self.multi_miller_loop(terms)
    }

    fn final_exponentiation(
        &mut self,
        f: &AssignedFq12<<G1Affine as halo2_proofs::arithmetic::CurveAffine>::Base, Fr>,
    ) -> AssignedFq12<<G1Affine as halo2_proofs::arithmetic::CurveAffine>::Base, Fr> {
        self.final_exponentiation(f)
    }
}
