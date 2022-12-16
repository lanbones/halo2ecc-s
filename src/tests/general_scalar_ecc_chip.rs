use crate::assign::AssignedPoint;
use crate::circuit::ecc_chip::{EccChipBaseOps, EccChipScalarOps};
use crate::circuit::integer_chip::IntegerChipOps;
use crate::context::Context;
use crate::context::GeneralScalarEccContext;
use crate::tests::{random_bls12_381_fr, run_circuit_on_bn256};
use crate::utils::field_to_bn;
use ark_std::{end_timer, start_timer};
use halo2_proofs::pairing::bls12_381::{G1Affine, G1};
use halo2_proofs::pairing::bn256::Fr;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_bls12_381_ecc_chip_over_bn256_fr() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let mut ctx = GeneralScalarEccContext::<G1Affine, Fr>::new(ctx);

    let mut points = vec![];
    let mut scalars = vec![];
    let mut acc = G1::identity();

    for _ in 0..50 {
        let a = random_bls12_381_fr();
        let b = random_bls12_381_fr();
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
        .map(|x| ctx.scalar_integer_ctx.assign_w(&field_to_bn(&x)))
        .collect::<Vec<_>>();
    let res: AssignedPoint<_, _> = ctx.msm(&points, &scalars);
    let res_expect: AssignedPoint<G1Affine, Fr> = ctx.assign_point(&acc);
    ctx.ecc_assert_equal(&res, &res_expect);

    end_timer!(timer);

    run_circuit_on_bn256(ctx.into(), 22);
}
