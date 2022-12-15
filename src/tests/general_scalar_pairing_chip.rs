use crate::assign::{AssignedCondition, AssignedG2Affine};
use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::circuit::fq12::{Fq12ChipOps, Fq2ChipOps};
use crate::circuit::pairing_chip::PairingChipOps;
use crate::context::{Context, GeneralScalarEccContext};
use crate::tests::run_circuit_on_bn256;
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

        let ab = pairing(&a, &b);

        let bx = ctx.fq2_assign_constant((b.x.c0, b.x.c1));
        let by = ctx.fq2_assign_constant((b.y.c0, b.y.c1));
        let b = AssignedG2Affine::new(
            bx,
            by,
            AssignedCondition(ctx.native_ctx.borrow_mut().assign_constant(Fr::zero())),
        );

        let ab0 = ctx.fq12_assign_constant((
            (
                (ab.0.c0.c0.c0, ab.0.c0.c0.c1),
                (ab.0.c0.c1.c0, ab.0.c0.c1.c1),
                (ab.0.c0.c2.c0, ab.0.c0.c2.c1),
            ),
            (
                (ab.0.c1.c0.c0, ab.0.c1.c0.c1),
                (ab.0.c1.c1.c0, ab.0.c1.c1.c1),
                (ab.0.c1.c2.c0, ab.0.c1.c2.c1),
            ),
        ));

        let a = ctx.assign_point(&a.to_curve());

        let ab1 = ctx.pairing(&[(&a, &b)]);

        ctx.fq12_assert_eq(&ab0, &ab1);

        run_circuit_on_bn256(ctx.into(), 22);
    }

    {
        let ctx = Rc::new(RefCell::new(Context::new()));
        let mut ctx = GeneralScalarEccContext::<G1Affine, Fr>::new(ctx);

        let a = G1::random(&mut OsRng);
        let b = G2Affine::from(G2::random(&mut OsRng));

        let bx = ctx.fq2_assign_constant((b.x.c0, b.x.c1));
        let by = ctx.fq2_assign_constant((b.y.c0, b.y.c1));
        let b = AssignedG2Affine::new(
            bx,
            by,
            AssignedCondition(ctx.native_ctx.borrow_mut().assign_constant(Fr::zero())),
        );

        let neg_a = ctx.assign_point(&-a);
        let a = ctx.assign_point(&a);

        ctx.check_pairing(&[(&a, &b), (&neg_a, &b)]);

        run_circuit_on_bn256(ctx.into(), 22);
    }
}
