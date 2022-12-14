use crate::circuit::base_chip::{BaseChip, BaseChipConfig};
use crate::circuit::range_chip::RangeChipConfig;
use crate::circuit::range_chip::{RangeChip, RangeChipOps};
use crate::context::Context;
use crate::context::Records;
use crate::range_info::{RangeInfo, MAX_CHUNKS, LIMBS, LIMB_BITS};
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
use std::marker::PhantomData;
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

        range_chip.init_table(&mut layouter, &RangeInfo::<N>::new::<W>())?;

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
    let ctx = Context::<Fq, Fr>::new_with_range_info();

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
        *ctx.range_offset = start;
        let common = &a & ((BigUint::from(1u64) << LIMB_BITS) - 1u64);
        let a_leading = &a >> ((LIMBS - 1) * LIMB_BITS);
        let b_leading = &b >> ((LIMBS - 1) * LIMB_BITS);
        let d_leading = &d >> ((LIMBS - 1) * LIMB_BITS);
        let r_leading = &r >> ((LIMBS - 1) * LIMB_BITS);
        let t = std::thread::spawn(move || {
            for _ in 0..step / (1 + MAX_CHUNKS) / 5 {
                ctx.assign_nonleading_limb(&common);
                ctx.assign_w_ceil_leading_limb(&a_leading);
                ctx.assign_w_ceil_leading_limb(&b_leading);
                ctx.assign_w_ceil_leading_limb(&r_leading);
                ctx.assign_d_leading_limb(&d_leading);
            }
            assert!(*ctx.range_offset <= start + step);
        });
        threads.push(t);
    }

    for t in threads {
        t.join().unwrap();
    }

    end_timer!(timer);

    const K: u32 = 19;
    let circuit = TestCircuit::<Fq, Fr> {
        records: Arc::try_unwrap(ctx.records).unwrap().into_inner().unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}
