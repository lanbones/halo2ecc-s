use crate::assign::AssignedPoint;
use crate::circuit::base_chip::{BaseChip, BaseChipConfig, BaseChipOps};
use crate::circuit::native_scalar_ecc_chip::EccChipOps;
use crate::circuit::range_chip::RangeChip;
use crate::circuit::range_chip::RangeChipConfig;
use crate::context::Records;
use crate::context::{Context, NativeScalarEccContext};
use crate::tests::random_fr;
use ark_std::{end_timer, start_timer};
use halo2_proofs::arithmetic::FieldExt;
use halo2_proofs::pairing::bn256::{Fq, Fr, G1Affine, G1};
use halo2_proofs::pairing::group::Group;
use halo2_proofs::{
    circuit::{Layouter, SimpleFloorPlanner},
    dev::MockProver,
    plonk::{Circuit, ConstraintSystem, Error},
};
use std::marker::PhantomData;
use std::sync::Arc;

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
fn test_native_ecc_chip() {
    let mut ctx = NativeScalarEccContext::<G1Affine>(Context::new_with_range_info());

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
        .map(|x| EccChipOps::<G1Affine>::assign_point(&mut ctx, x))
        .collect::<Vec<_>>();
    let scalars = scalars
        .into_iter()
        .map(|x| ctx.0.assign(x))
        .collect::<Vec<_>>();
    let res: AssignedPoint<_, _> = ctx.msm(&points, &scalars);
    let res_expect: AssignedPoint<G1Affine, Fr> = ctx.assign_point(&acc);
    ctx.ecc_assert_equal(&res, &res_expect);

    end_timer!(timer);

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
