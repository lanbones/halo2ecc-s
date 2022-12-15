use crate::assign::{AssignedCondition, AssignedG2Affine};
use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::circuit::fq12::Fq2ChipOps;
use crate::circuit::pairing_chip::PairingChipOps;
use crate::context::IntegerContext;
use crate::context::{Context, NativeScalarEccContext};
use crate::tests::run_circuit_on_bn256;
use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::pairing::bn256::{Fr, G1Affine, G2Affine, G1, G2};
use halo2_proofs::pairing::group::Group;
use rand::thread_rng;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_native_pairing_chip() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
    let mut ctx = NativeScalarEccContext::<G1Affine>(ctx);

    let a = G1::random(&mut thread_rng());
    let b = G2Affine::from(G2::random(&mut thread_rng()));

    let bx = ctx.fq2_assign_constant(
        b.coordinates().unwrap().x().c0,
        b.coordinates().unwrap().x().c1,
    );
    let by = ctx.fq2_assign_constant(
        b.coordinates().unwrap().y().c0,
        b.coordinates().unwrap().y().c1,
    );
    let b = AssignedG2Affine::new(
        bx,
        by,
        AssignedCondition(ctx.0.ctx.borrow_mut().assign_constant(Fr::zero())),
    );
    let neg_a = ctx.assign_point(&-a);
    let a = ctx.assign_point(&a);

    ctx.check_pairing(&[(&a, &b), (&neg_a, &b)]);

    run_circuit_on_bn256(Rc::try_unwrap(ctx.0.ctx).unwrap().into_inner(), 22);
}
