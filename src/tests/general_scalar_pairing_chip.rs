use crate::assign::{AssignedCondition, AssignedG2Affine};
use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::circuit::fq12::{Fq12ChipOps, Fq2ChipOps};
use crate::circuit::pairing_chip::PairingChipOps;
use crate::context::{Context, GeneralScalarEccContext};
use crate::tests::run_circuit_on_bn256;
use halo2_proofs::arithmetic::Field;
use halo2_proofs::pairing::bls12_381::pairing;
use halo2_proofs::pairing::bls12_381::{G1Affine, G2Affine, G1, G2};
use halo2_proofs::pairing::bn256::Fr;
use halo2_proofs::pairing::group::prime::PrimeCurveAffine;
use halo2_proofs::pairing::group::Group;
use rand::rngs::OsRng;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_bls12_381_pairing_chip_over_bn256_fr() {
    {
        let ctx = Rc::new(RefCell::new(Context::new()));
        let mut ctx = GeneralScalarEccContext::<G1Affine, Fr>::new(ctx);

        let a = G1::random(&mut OsRng).into();
        let b = G2Affine::from(G2::random(&mut OsRng));
        let c = G1::random(&mut OsRng).into();
        let d = G2Affine::from(G2::random(&mut OsRng));

        let ab = pairing(&a, &b);
        let cd = pairing(&c, &d);

        let abcd = ab + cd;

        let bx = ctx.fq2_assign_constant((b.x.c0, b.x.c1));
        let by = ctx.fq2_assign_constant((b.y.c0, b.y.c1));
        let b = AssignedG2Affine::new(
            bx,
            by,
            AssignedCondition(ctx.native_ctx.borrow_mut().assign_constant(Fr::zero())),
        );

        let dx = ctx.fq2_assign_constant((d.x.c0, d.x.c1));
        let dy = ctx.fq2_assign_constant((d.y.c0, d.y.c1));
        let d = AssignedG2Affine::new(
            dx,
            dy,
            AssignedCondition(ctx.native_ctx.borrow_mut().assign_constant(Fr::zero())),
        );

        let abcd0 = ctx.fq12_assign_constant((
            (
                (abcd.0.c0.c0.c0, abcd.0.c0.c0.c1),
                (abcd.0.c0.c1.c0, abcd.0.c0.c1.c1),
                (abcd.0.c0.c2.c0, abcd.0.c0.c2.c1),
            ),
            (
                (abcd.0.c1.c0.c0, abcd.0.c1.c0.c1),
                (abcd.0.c1.c1.c0, abcd.0.c1.c1.c1),
                (abcd.0.c1.c2.c0, abcd.0.c1.c2.c1),
            ),
        ));

        let a = ctx.assign_point(&a.to_curve());
        let c = ctx.assign_point(&c.to_curve());

        let abcd1 = ctx.pairing(&[(&a, &b), (&c, &d)]);

        ctx.fq12_assert_eq(&abcd0, &abcd1);

        run_circuit_on_bn256(ctx.into(), 22);
    }

    {
        let ctx = Rc::new(RefCell::new(Context::new()));
        let mut ctx = GeneralScalarEccContext::<G1Affine, Fr>::new(ctx);

        let a = G1::random(&mut OsRng);
        let b = G2Affine::from(G2::random(&mut OsRng));
        let c = halo2_proofs::pairing::bls12_381::Fr::random(&mut OsRng);
        let ac = a * c;
        let bc = G2Affine::from(b * c);

        let bx = ctx.fq2_assign_constant((b.x.c0, b.x.c1));
        let by = ctx.fq2_assign_constant((b.y.c0, b.y.c1));
        let b = AssignedG2Affine::new(
            bx,
            by,
            AssignedCondition(ctx.native_ctx.borrow_mut().assign_constant(Fr::zero())),
        );

        let bcx = ctx.fq2_assign_constant((bc.x.c0, bc.x.c1));
        let bcy = ctx.fq2_assign_constant((bc.y.c0, bc.y.c1));
        let bc = AssignedG2Affine::new(
            bcx,
            bcy,
            AssignedCondition(ctx.native_ctx.borrow_mut().assign_constant(Fr::zero())),
        );

        let neg_a = ctx.assign_point(&-a);
        let ac = ctx.assign_point(&ac);

        ctx.check_pairing(&[(&ac, &b), (&neg_a, &bc)]);

        run_circuit_on_bn256(ctx.into(), 22);
    }
}
