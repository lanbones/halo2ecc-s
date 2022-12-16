use std::sync::Arc;

use ark_std::{end_timer, start_timer};
use halo2_proofs::{
    arithmetic::{BaseExt, FieldExt},
    circuit::{Layouter, SimpleFloorPlanner},
    dev::MockProver,
    pairing::bn256::Fr,
    plonk::{Circuit, ConstraintSystem, Error},
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use crate::{
    circuit::{
        base_chip::{BaseChip, BaseChipConfig},
        range_chip::{RangeChip, RangeChipConfig},
    },
    context::{Context, Records},
};

pub mod base_chip;
pub mod general_scalar_ecc_chip;
pub mod general_scalar_pairing_chip;
pub mod integer_chip;
pub mod native_scalar_ecc_chip;
pub mod native_scalar_pairing_chip;
pub mod range_chip;

fn random<N: BaseExt>() -> N {
    let seed = chrono::offset::Utc::now()
        .timestamp_nanos()
        .try_into()
        .unwrap();
    let rng = XorShiftRng::seed_from_u64(seed);
    N::random(rng)
}

fn random_fr() -> halo2_proofs::pairing::bn256::Fr {
    random::<halo2_proofs::pairing::bn256::Fr>()
}

fn random_fq() -> halo2_proofs::pairing::bn256::Fq {
    random::<halo2_proofs::pairing::bn256::Fq>()
}

fn random_bls12_381_fr() -> halo2_proofs::pairing::bls12_381::Fr {
    random::<halo2_proofs::pairing::bls12_381::Fr>()
}

fn random_bls12_381_fq() -> halo2_proofs::pairing::bls12_381::Fq {
    random::<halo2_proofs::pairing::bls12_381::Fq>()
}

#[derive(Clone)]
struct TestChipConfig {
    base_chip_config: BaseChipConfig,
    range_chip_config: RangeChipConfig,
}

#[derive(Default, Clone)]
struct TestCircuit<N: FieldExt> {
    records: Records<N>,
}

impl<N: FieldExt> Circuit<N> for TestCircuit<N> {
    type Config = TestChipConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<N>) -> Self::Config {
        let base_chip_config = BaseChip::configure(meta);
        let range_chip_config = RangeChip::<N>::configure(meta);
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
        let range_chip = RangeChip::<N>::new(config.range_chip_config);

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

pub fn run_circuit_on_bn256(ctx: Context<Fr>, k: u32) {
    println!("offset {} {}", ctx.range_offset, ctx.base_offset);

    let circuit = TestCircuit::<Fr> {
        records: Arc::try_unwrap(ctx.records).unwrap().into_inner().unwrap(),
    };

    let prover = match MockProver::run(k, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}
