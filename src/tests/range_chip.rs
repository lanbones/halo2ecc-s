use crate::circuit::range_chip::RangeChipOps;
use crate::circuit::range_chip::MAX_CHUNKS;
use crate::context::{Context, IntegerContext};
use crate::tests::random_fq;
use crate::tests::run_circuit_on_bn256;
use crate::utils::field_to_bn;
use ark_std::{end_timer, start_timer};
use halo2_proofs::pairing::bn256::{Fq, Fr};
use num_bigint::BigUint;
use num_integer::Integer;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_range_chip() {
    let ctx = Context::new();

    let a = field_to_bn(&random_fq());
    let b = field_to_bn(&random_fq());
    let fq_modulus = field_to_bn(&-Fq::one()) + 1u64;
    let (d, r) = (&a * &b).div_rem(&fq_modulus);

    let timer = start_timer!(|| "setup");
    let c = (1 + MAX_CHUNKS) * 5 * 330;
    let n = 10;

    let mut threads = vec![];
    for i in 0..n {
        let step = c / n;
        let start = i * step * 2;
        let mut ctx = ctx.clone();
        ctx.range_offset = start as usize;
        let a = a.clone();
        let b = b.clone();
        let d = d.clone();
        let r = r.clone();
        let t = std::thread::spawn(move || {
            let mut ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(Rc::new(
                RefCell::new(ctx),
            ));

            let limbs = ctx.info().limbs;
            let non_leading_bits = (limbs - 1) * ctx.info().limb_bits;
            let common = &a & ((BigUint::from(1u64) << ctx.info().limbs) - 1u64);
            let a_leading = &a >> non_leading_bits;
            let b_leading = &b >> non_leading_bits;
            let d_leading = &d >> non_leading_bits;
            let r_leading = &r >> non_leading_bits;
            for _ in 0..step / (1 + MAX_CHUNKS) / 5 {
                ctx.assign_nonleading_limb(&common);
                ctx.assign_w_ceil_leading_limb(&a_leading);
                ctx.assign_w_ceil_leading_limb(&b_leading);
                ctx.assign_w_ceil_leading_limb(&r_leading);
                ctx.assign_d_leading_limb(&d_leading);
            }
            assert!(ctx.ctx.borrow().range_offset as u64 <= start + step);
        });
        threads.push(t);
    }

    for t in threads {
        t.join().unwrap();
    }

    end_timer!(timer);

    run_circuit_on_bn256(ctx, 20);
}
