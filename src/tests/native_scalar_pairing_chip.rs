use crate::assign::{AssignedCondition, AssignedG2Affine, AssignedPoint};
use crate::circuit::base_chip::{BaseChip, BaseChipConfig, BaseChipOps};
use crate::circuit::native_scalar_ecc_chip::EccChipOps;
use crate::circuit::range_chip::RangeChip;
use crate::circuit::range_chip::RangeChipConfig;
use crate::context::Records;
use crate::context::{Context, NativeScalarEccContext};
use crate::range_info::RangeInfo;
use crate::tests::random_fr;
use ark_std::{end_timer, start_timer};
use halo2_proofs::arithmetic::{CurveAffine, FieldExt};
use halo2_proofs::pairing::bn256::{Bn256, Fq, Fr, G1Affine, G2Affine, G1, G2};
use halo2_proofs::pairing::group::Group;
use halo2_proofs::plonk::{create_proof, keygen_pk, keygen_vk};
use halo2_proofs::poly::commitment::Params;
use halo2_proofs::transcript::{Blake2bRead, Blake2bWrite, Challenge255};
use halo2_proofs::{
    circuit::{Layouter, SimpleFloorPlanner},
    dev::MockProver,
    plonk::{Circuit, ConstraintSystem, Error},
};
use rand::rngs::OsRng;
use rand::thread_rng;
use rand_xorshift::XorShiftRng;
use std::marker::PhantomData;
use std::sync::Arc;
use std::thread;

#[derive(Clone)]
struct TestChipConfig {
    base_chip_config: BaseChipConfig,
    range_chip_config: RangeChipConfig,
}

#[derive(Default, Clone)]
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
fn test_native_pairing_chip() {
    let mut ctx = NativeScalarEccContext::<G1Affine>(Context::new_with_range_info());

    let a = G1::random(&mut thread_rng());
    let b = G2Affine::from(G2::random(&mut thread_rng()));

    let bx = ctx.0.fq2_assign_constant(
        b.coordinates().unwrap().x().c0,
        b.coordinates().unwrap().x().c1,
    );
    let by = ctx.0.fq2_assign_constant(
        b.coordinates().unwrap().y().c0,
        b.coordinates().unwrap().y().c1,
    );
    let b = AssignedG2Affine::new(bx, by, AssignedCondition(ctx.0.assign_constant(Fr::zero())));
    let neg_a = ctx.assign_point(&-a);
    let a = ctx.assign_point(&a);

    ctx.check_pairing(&[(&a, &b), (&neg_a, &b)]);

    println!("offset {} {}", ctx.0.range_offset, ctx.0.base_offset);

    const K: u32 = 22;
    let circuit = TestCircuit::<Fq, Fr> {
        records: Arc::try_unwrap(ctx.0.records)
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
