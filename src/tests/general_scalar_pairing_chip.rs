use crate::assign::{AssignedCondition, AssignedG2Affine};
use crate::context::{Context, GeneralScalarEccContext};
use crate::tests::run_circuit_on_bn256;
use rand::thread_rng;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_bls12_381_on_bn256_pairing_chip() {
    use crate::circuit::pairing_chip::PairingChipOps;
    use halo2_proofs::pairing::bls12_381::G1Affine;
    use halo2_proofs::pairing::bls12_381::G2Affine;
    use halo2_proofs::pairing::bls12_381::G1;
    use halo2_proofs::pairing::bls12_381::G2;
    use halo2_proofs::pairing::bn256::Fr;
    use halo2_proofs::pairing::group::Group;
    use crate::circuit::fq12::Fq2ChipOps;
    use crate::circuit::ecc_chip::EccChipBaseOps;
    use crate::circuit::ecc_chip::EccBaseIntegerChipWrapper;

    let ctx = Rc::new(RefCell::new(Context::new()));
    let mut ctx = GeneralScalarEccContext::<G1Affine, Fr>::new(ctx);

    let a = G1::random(&mut thread_rng());
    let b = G2Affine::from(G2::random(&mut thread_rng()));

    let bx = ctx.fq2_assign_constant(b.x.c0, b.x.c1);
    let by = ctx.fq2_assign_constant(b.y.c0, b.y.c1);
    let b = AssignedG2Affine::new(
        bx,
        by,
        AssignedCondition(ctx.base_integer_chip().base_chip().assign_constant(Fr::zero())),
    );
    let neg_a = ctx.assign_point(&-a);
    let a = ctx.assign_point(&a);

    ctx.check_pairing(&[(&a, &b), (&neg_a, &b)]);

    run_circuit_on_bn256(ctx.into(), 22);
}
