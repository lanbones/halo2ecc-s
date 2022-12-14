use crate::circuit::base_chip::{BaseChip, BaseChipConfig, BaseChipOps, MUL_COLUMNS, VAR_COLUMNS};
use crate::context::Context;
use crate::context::Records;
use crate::pair;
use crate::tests::random_fr;
use ark_std::{end_timer, start_timer};
use halo2_proofs::arithmetic::FieldExt;
use halo2_proofs::pairing::bn256::Fr;
use halo2_proofs::{
    circuit::{Layouter, SimpleFloorPlanner},
    dev::MockProver,
    plonk::{Circuit, ConstraintSystem, Error},
};
use std::marker::PhantomData;
use std::sync::Arc;

#[derive(Clone)]
struct TestBaseChipConfig {
    base_chip_config: BaseChipConfig,
}

#[derive(Default)]
struct TestCircuit<N: FieldExt> {
    records: Records<N>,
    _phantom: PhantomData<N>,
}

impl<N: FieldExt> Circuit<N> for TestCircuit<N> {
    type Config = TestBaseChipConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<N>) -> Self::Config {
        let base_chip_config = BaseChip::configure(meta);
        TestBaseChipConfig { base_chip_config }
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<N>,
    ) -> Result<(), Error> {
        let base_chip = BaseChip::new(config.base_chip_config);

        layouter.assign_region(
            || "base",
            |mut region| {
                let timer = start_timer!(|| "assign");
                self.records.assign_all_in_base(&mut region, &base_chip)?;
                end_timer!(timer);
                Ok(())
            },
        )?;

        Ok(())
    }
}

#[test]
fn test_one_line_st() {
    let vars = [(); VAR_COLUMNS].map(|_| random_fr());
    let coeffs = [(); VAR_COLUMNS].map(|_| random_fr());
    let muls_coeffs = [(); MUL_COLUMNS].map(|_| random_fr());
    let next_var = random_fr();
    let next_coeff = random_fr();

    let result: Fr = {
        let mut result = Fr::zero();
        for i in 0..VAR_COLUMNS {
            result = result + vars[i] * coeffs[i]
        }
        for i in 0..MUL_COLUMNS {
            result = result + muls_coeffs[i] * vars[i * 2] * vars[i * 2 + 1]
        }
        result + next_var * next_coeff
    };

    let ctx = Context::<Fr>::new();
    {
        let ctx = &mut ctx.clone();

        let timer = start_timer!(|| "setup");
        for _ in 0..10000 {
            ctx.one_line(
                (0..VAR_COLUMNS)
                    .map(|i| pair!(vars[i], coeffs[i]))
                    .collect(),
                Some(-result),
                (muls_coeffs.to_vec(), Some(next_coeff)),
            );

            ctx.one_line_with_last(vec![], pair!(next_var, Fr::zero()), None, (vec![], None));
        }
        end_timer!(timer);
    }

    const K: u32 = 16;
    let circuit = TestCircuit::<Fr> {
        records: Arc::try_unwrap(ctx.records).unwrap().into_inner().unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}

#[test]
fn test_one_line_mt() {
    let vars = [(); VAR_COLUMNS].map(|_| random_fr());
    let coeffs = [(); VAR_COLUMNS].map(|_| random_fr());
    let muls_coeffs = [(); MUL_COLUMNS].map(|_| random_fr());
    let next_var = random_fr();
    let next_coeff = random_fr();

    let result = {
        let mut result = Fr::zero();
        for i in 0..VAR_COLUMNS {
            result = result + vars[i] * coeffs[i]
        }
        for i in 0..MUL_COLUMNS {
            result = result + muls_coeffs[i] * vars[i * 2] * vars[i * 2 + 1]
        }
        result + next_var * next_coeff
    };

    let mut threads = vec![];

    let timer = start_timer!(|| "setup");
    let c = 10000;
    let n = 10;
    let ctx = Context::<Fr>::new();
    for i in 0..n {
        let step = c / n;
        let start = i * step * 2;
        let mut ctx = ctx.clone();
        *ctx.base_offset = start;
        let t = std::thread::spawn(move || {
            for _ in 0..step {
                ctx.one_line(
                    (0..VAR_COLUMNS)
                        .map(|i| pair!(vars[i], coeffs[i]))
                        .collect(),
                    Some(-result),
                    (muls_coeffs.to_vec(), Some(next_coeff)),
                );

                ctx.one_line_with_last(vec![], pair!(next_var, Fr::zero()), None, (vec![], None));
            }
        });
        threads.push(t);
    }

    for t in threads {
        t.join().unwrap();
    }
    end_timer!(timer);

    const K: u32 = 16;
    let circuit = TestCircuit::<Fr> {
        records: Arc::try_unwrap(ctx.records).unwrap().into_inner().unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}
