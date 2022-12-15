/*
  The implementation is ported from https://github.com/zkcrypto/bls12_381
*/

use halo2_proofs::{
    arithmetic::{CurveAffine, FieldExt},
    pairing::{
        bls12_381::{Fq, G1Affine},
        bn256::Fr,
    },
};

use crate::{
    assign::{
        AssignedFq12, AssignedFq2, AssignedFq6, AssignedG1Affine, AssignedG2Affine,
        AssignedG2Prepared,
    },
    context::GeneralScalarEccContext,
};

use super::{
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
        let c0 = self.base_integer_chip().int_sub(&a.0, &a.1);
        let c1 = self.base_integer_chip().int_add(&a.0, &a.1);

        (c0, c1)
    }

    fn fq2_frobenius_map(&mut self, x: &AssignedFq2<Fq, N>, _power: usize) -> AssignedFq2<Fq, N> {
        self.fq2_conjugate(x)
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

        let coeff_c1 = self.fq2_assign_constant((
            Fq::zero(),
            Fq::from_raw_unchecked([
                0xcd03_c9e4_8671_f071,
                0x5dab_2246_1fcd_a5d2,
                0x5870_42af_d385_1b95,
                0x8eb6_0ebe_01ba_cb9e,
                0x03f9_7d6e_83d0_50d2,
                0x18f0_2065_5463_8741,
            ]),
        ));
        let c1 = self.fq2_mul(&c1, &coeff_c1);
        let coeff_c2 = self.fq2_assign_constant((
            Fq::zero(),
            Fq::from_raw_unchecked([
                0x890d_c9e4_8675_45c3,
                0x2af3_2253_3285_a5d5,
                0x5088_0866_309b_7e2c,
                0xa20d_1b8c_7e88_1024,
                0x14e4_f04f_e2db_9068,
                0x14e5_6d3f_1564_853a,
            ]),
        ));
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

        let coeff = self.fq2_assign_constant((
            Fq::from_raw_unchecked([
                0x0708_9552_b319_d465,
                0xc669_5f92_b50a_8313,
                0x97e8_3ccc_d117_228f,
                0xa35b_aeca_b2dc_29ee,
                0x1ce3_93ea_5daa_ce4d,
                0x08f2_220f_b0fb_66eb,
            ]),
            Fq::from_raw_unchecked([
                0xb2f6_6aad_4ce5_d646,
                0x5842_a06b_fc49_7cec,
                0xcf48_95d4_2599_d394,
                0xc11b_9cba_40a8_e8d0,
                0x2e38_13cb_e5a0_de89,
                0x110e_efda_8884_7faf,
            ]),
        ));
        let c1c0 = self.fq2_mul(&c1.0, &coeff);
        let c1c1 = self.fq2_mul(&c1.1, &coeff);
        let c1c2 = self.fq2_mul(&c1.2, &coeff);

        (c0, (c1c0, c1c1, c1c2))
    }
}

impl Fq2ChipOps<Fq, Fr> for GeneralScalarEccContext<G1Affine, Fr> {}
impl Fq6ChipOps<Fq, Fr> for GeneralScalarEccContext<G1Affine, Fr> {}
impl Fq12ChipOps<Fq, Fr> for GeneralScalarEccContext<G1Affine, Fr> {}

impl GeneralScalarEccContext<G1Affine, Fr> {
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

    fn cyclotomic_square(&mut self, _f: &AssignedFq12<Fq, Fr>) -> AssignedFq12<Fq, Fr> {
        unimplemented!()
    }

    fn cycolotomic_exp(&mut self, _f: &AssignedFq12<Fq, Fr>) -> AssignedFq12<Fq, Fr> {
        unimplemented!()
    }
}

const BLS_X: u64 = 0xd201_0000_0001_0000;

impl PairingChipOps<G1Affine, Fr> for GeneralScalarEccContext<G1Affine, Fr> {
    fn prepare_g2(
        &mut self,
        g2: &AssignedG2Affine<G1Affine, Fr>,
    ) -> AssignedG2Prepared<G1Affine, Fr> {
        let mut f = self.g2affine_to_g2(g2);
        let mut coeffs = vec![];

        let mut found_one = false;
        for i in (0..64).rev().map(|b| (((BLS_X >> 1) >> b) & 1) == 1) {
            if !found_one {
                found_one = i;
                continue;
            }

            coeffs.push(self.doubling_step(&mut f));

            if i {
                coeffs.push(self.addition_step(&mut f, g2));
            }
        }

        coeffs.push(self.doubling_step(&mut f));

        AssignedG2Prepared::new(coeffs)
    }

    fn multi_miller_loop(
        &mut self,
        terms: &[(
            &AssignedG1Affine<G1Affine, Fr>,
            &AssignedG2Prepared<G1Affine, Fr>,
        )],
    ) -> AssignedFq12<<G1Affine as halo2_proofs::arithmetic::CurveAffine>::Base, Fr> {
        let mut pairs = vec![];
        for &(p, q) in terms {
            // not support identity
            self.base_integer_chip().base_chip().assert_false(&p.z);
            pairs.push((p, q.coeffs.iter()));
        }

        let mut f = self.fq12_assign_one();

        let mut found_one = false;
        for i in (0..64).rev().map(|b| (((BLS_X >> 1) >> b) & 1) == 1) {
            if !found_one {
                found_one = i;
                continue;
            }

            for &mut (p, ref mut coeffs) in &mut pairs {
                f = self.ell(&f, coeffs.next().unwrap(), &p);
            }

            if i {
                for &mut (p, ref mut coeffs) in &mut pairs {
                    f = self.ell(&f, coeffs.next().unwrap(), &p);
                }
            }

            f = self.fq12_square(&f);
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            f = self.ell(&f, coeffs.next().unwrap(), &p);
        }

        f = self.fq12_conjugate(&f);

        f
    }

    fn final_exponentiation(
        &mut self,
        f: &AssignedFq12<<G1Affine as CurveAffine>::Base, Fr>,
    ) -> AssignedFq12<<G1Affine as CurveAffine>::Base, Fr> {
        const POWER_HOLDER: usize = 1;

        let mut t0 = self.fq12_frobenius_map(f, POWER_HOLDER);
        t0 = self.fq12_frobenius_map(&t0, POWER_HOLDER);
        t0 = self.fq12_frobenius_map(&t0, POWER_HOLDER);
        t0 = self.fq12_frobenius_map(&t0, POWER_HOLDER);
        t0 = self.fq12_frobenius_map(&t0, POWER_HOLDER);
        t0 = self.fq12_frobenius_map(&t0, POWER_HOLDER);

        let mut t1 = self.fq12_unsafe_invert(&f);

        let mut t2 = self.fq12_mul(&t0, &t1);
        t1 = t2.clone();

        t2 = self.fq12_frobenius_map(&t2, POWER_HOLDER);
        t2 = self.fq12_frobenius_map(&t2, POWER_HOLDER);

        t2 = self.fq12_mul(&t2, &t1);
        t1 = self.cyclotomic_square(&t2);
        t1 = self.fq12_conjugate(&t1);
        let mut t3 = self.cycolotomic_exp(&t2);
        let mut t4 = self.cyclotomic_square(&t3);
        let mut t5 = self.fq12_mul(&t1, &t3);
        t1 = self.cycolotomic_exp(&t5);
        t0 = self.cycolotomic_exp(&t1);
        let mut t6 = self.cycolotomic_exp(&t0);
        t6 = self.fq12_mul(&t6, &t4);
        t4 = self.cycolotomic_exp(&t6);
        t5 = self.fq12_conjugate(&t5);
        let t = self.fq12_mul(&t5, &t2);
        t4 = self.fq12_mul(&t4, &t);
        t5 = self.fq12_conjugate(&t2);
        t1 = self.fq12_mul(&t1, &t2);
        for _ in 0..3 {
            t1 = self.fq12_frobenius_map(&t1, POWER_HOLDER)
        }
        t6 = self.fq12_mul(&t6, &t5);
        t6 = self.fq12_frobenius_map(&t6, POWER_HOLDER);
        t3 = self.fq12_mul(&t3, &t0);
        for _ in 0..2 {
            t3 = self.fq12_frobenius_map(&t3, POWER_HOLDER)
        }
        t3 = self.fq12_mul(&t3, &t1);
        t3 = self.fq12_mul(&t3, &t6);

        self.fq12_mul(&t3, &t4)
    }
}
