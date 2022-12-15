use crate::assign::AssignedPoint;
use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::ecc_chip::{EccChipBaseOps, EccChipScalarOps};
use crate::context::IntegerContext;
use crate::context::{Context, NativeScalarEccContext};
use crate::tests::{random_fr, run_circuit_on_bn256};
use ark_std::{end_timer, start_timer};
use halo2_proofs::pairing::bn256::{Fr, G1Affine, G1};
use halo2_proofs::pairing::group::Group;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_native_ecc_chip() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
    let mut ctx = NativeScalarEccContext(ctx);

    let mut points = vec![];
    let mut scalars = vec![];
    let mut acc = G1::identity();

    for _ in 0..100 {
        let a = random_fr();
        let b = random_fr();
        let p = G1::generator() * a;
        acc = acc + p * b;
        points.push(p);
        scalars.push(b);
    }

    let timer = start_timer!(|| "setup");

    let points = points
        .iter()
        .map(|x| ctx.assign_point(x))
        .collect::<Vec<_>>();
    let scalars = scalars
        .into_iter()
        .map(|x| ctx.0.ctx.borrow_mut().assign(x))
        .collect::<Vec<_>>();
    let res: AssignedPoint<_, _> = ctx.msm(&points, &scalars);
    let res_expect: AssignedPoint<G1Affine, Fr> = ctx.assign_point(&acc);
    ctx.ecc_assert_equal(&res, &res_expect);

    end_timer!(timer);

    run_circuit_on_bn256(Rc::try_unwrap(ctx.0.ctx).unwrap().into_inner(), 22);
}
