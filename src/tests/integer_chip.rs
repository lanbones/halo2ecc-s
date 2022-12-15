use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::integer_chip::IntegerChipOps;
use crate::context::{Context, IntegerContext};
use crate::tests::{random_bls12_381_fq, random_bls12_381_fr, random_fq, run_circuit_on_bn256};
use crate::utils::field_to_bn;
use halo2_proofs::arithmetic::Field;
use halo2_proofs::pairing::bn256::Fr;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_integer_chip_st() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let mut ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);

    let a = random_fq();
    let mut b = random_fq();

    let mut timeout = 3;
    while b.is_zero().into() && timeout > 0 {
        b = random_fq();
        timeout -= 1;
    }
    assert!(timeout > 0);

    let c = a + b;
    let d = a - b;
    let e = a * b;
    let f = a * b.invert().unwrap();

    let a = ctx.assign_w(&field_to_bn(&a));
    let b = ctx.assign_w(&field_to_bn(&b));

    let c1 = ctx.assign_w(&field_to_bn(&c));
    let c2 = ctx.int_add(&a, &b);
    ctx.assert_int_equal(&c1, &c2);

    let d1 = ctx.assign_w(&field_to_bn(&d));
    let d2 = ctx.int_sub(&a, &b);
    ctx.assert_int_equal(&d1, &d2);

    let e1 = ctx.assign_w(&field_to_bn(&e));
    let e2 = ctx.int_mul(&a, &b);
    ctx.assert_int_equal(&e1, &e2);

    let f1 = ctx.assign_w(&field_to_bn(&f));
    let f2 = ctx.int_div(&a, &b).1;
    ctx.assert_int_equal(&f1, &f2);

    let zero = ctx.int_sub(&a, &a);
    let (g1, _) = ctx.int_div(&a, &zero);
    ctx.ctx.borrow_mut().assert_true(&g1);

    run_circuit_on_bn256(ctx.into(), 20);
}

#[test]
fn test_bls12_381_fq_integer_chip_st() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let mut ctx = IntegerContext::<halo2_proofs::pairing::bls12_381::Fq, Fr>::new(ctx);

    for _ in 0..1000 {
        let a = random_bls12_381_fq();
        let b = random_bls12_381_fq();
        let ab = a * b;

        let a = ctx.assign_w(&field_to_bn(&a));
        let b = ctx.assign_w(&field_to_bn(&b));
        let ab0 = ctx.assign_w(&field_to_bn(&ab));

        let ab1 = ctx.int_mul(&a, &b);

        ctx.assert_int_equal(&ab0, &ab1);
    }

    run_circuit_on_bn256(ctx.into(), 20);
}

#[test]
fn test_bls12_381_fr_integer_chip_st() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let mut ctx = IntegerContext::<halo2_proofs::pairing::bls12_381::Fr, Fr>::new(ctx);

    for _ in 0..1000 {
        let a = random_bls12_381_fr();
        let b = random_bls12_381_fr();
        let ab = a * b;

        let a = ctx.assign_w(&field_to_bn(&a));
        let b = ctx.assign_w(&field_to_bn(&b));
        let ab0 = ctx.assign_w(&field_to_bn(&ab));

        let ab1 = ctx.int_mul(&a, &b);

        ctx.assert_int_equal(&ab0, &ab1);
    }

    run_circuit_on_bn256(ctx.into(), 20);
}
