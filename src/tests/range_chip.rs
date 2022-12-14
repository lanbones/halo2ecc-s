use crate::circuit::base_chip::{BaseChip, BaseChipConfig};
use crate::circuit::range_chip::{RangeChip, RangeChipOps};
use crate::circuit::range_chip::{RangeChipConfig, MAX_CHUNKS};
use crate::context::Records;
use crate::context::{Context, IntegerContext};
use crate::tests::random_fq;
use crate::utils::field_to_bn;
use ark_std::{end_timer, start_timer};
use halo2_proofs::arithmetic::FieldExt;
use halo2_proofs::pairing::bn256::{Fq, Fr};
use halo2_proofs::{
    circuit::{Layouter, SimpleFloorPlanner},
    dev::MockProver,
    plonk::{Circuit, ConstraintSystem, Error},
};
use num_bigint::BigUint;
use num_integer::Integer;
use std::cell::RefCell;
use std::marker::PhantomData;
use std::rc::Rc;
use std::sync::Arc;

#[derive(Clone)]
struct TestChipConfig {
    base_chip_config: BaseChipConfig,
    range_chip_config: RangeChipConfig,
}

#[derive(Default)]
struct TestCircuit<W: FieldExt, N: FieldExt> {
    records: Records<N>,
    _phantom: PhantomData<W>,
}

impl<W: FieldExt, N: FieldExt> Circuit<N> for TestCircuit<W, N> {
    type Config = TestChipConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<N>) -> Self::Config {
        let base_chip_config = BaseChip::configure(meta);
        let range_chip_config = RangeChip::<W, N>::configure(meta);
        TestChipConfig {
            base_chip_config,
            range_chip_config,
        }
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<N>,
    ) -> Result<(), Error> {
        let base_chip = BaseChip::new(config.base_chip_config);
        let range_chip = RangeChip::<W, N>::new(config.range_chip_config);

        range_chip.init_table(&mut layouter)?;

        layouter.assign_region(
            || "base",
            |mut region| {
                let timer = start_timer!(|| "assign");
                self.records
                    .assign_all(&mut region, &base_chip, &range_chip)?;
                end_timer!(timer);
                Ok(())
            },
        )?;

        Ok(())
    }
}

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

    const K: u32 = 20;
    let circuit = TestCircuit::<Fq, Fr> {
        records: Arc::try_unwrap(ctx.records)
            .unwrap()
            .into_inner()
            .unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}
