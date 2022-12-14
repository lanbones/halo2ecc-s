use crate::assign::AssignedPoint;
use crate::circuit::base_chip::{BaseChip, BaseChipConfig, BaseChipOps};
use crate::circuit::native_ecc_chip::EccChipOps;
use crate::circuit::range_chip::RangeChip;
use crate::circuit::range_chip::RangeChipConfig;
use crate::context::EccContext;
use crate::context::Records;
use crate::range_info::RangeInfo;
use crate::tests::random_fr;
use ark_std::{end_timer, start_timer};
use halo2_proofs::arithmetic::FieldExt;
use halo2_proofs::pairing::bn256::{Bn256, Fq, Fr, G1Affine, G1};
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
fn test_native_ecc_chip() {
    let mut ctx = EccContext::<G1Affine>::new_with_range_info();

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
        .map(|x| ctx.assign(x))
        .collect::<Vec<_>>();
    let res: AssignedPoint<_, _> = ctx.msm(&points, &scalars);
    let res_expect: AssignedPoint<G1Affine, Fr> = ctx.assign_point(&acc);
    ctx.ecc_assert_equal(&res, &res_expect);

    end_timer!(timer);

    println!("offset {} {}", ctx.range_offset, ctx.base_offset);

    const K: u32 = 22;
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

#[test]
fn bench_native_ecc_chip() {
    let mut ctx = EccContext::<G1Affine>::new_with_range_info();

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
        .map(|x| ctx.assign(x))
        .collect::<Vec<_>>();
    let res: AssignedPoint<_, _> = ctx.msm(&points, &scalars);
    let res_expect: AssignedPoint<G1Affine, Fr> = ctx.assign_point(&acc);
    ctx.ecc_assert_equal(&res, &res_expect);

    end_timer!(timer);

    const K: u32 = 22;
    let circuit = TestCircuit::<Fq, Fr> {
        records: Arc::try_unwrap(ctx.records).unwrap().into_inner().unwrap(),
        _phantom: PhantomData,
    };

    let timer = start_timer!(|| "setup params");
    let params = Params::<G1Affine>::unsafe_setup::<Bn256>(K);
    end_timer!(timer);

    let timer = start_timer!(|| "gen vk");
    let vk = keygen_vk(&params, &circuit).expect("keygen_vk should not fail");
    end_timer!(timer);

    let timer = start_timer!(|| "gen pk");
    let pk = keygen_pk(&params, vk, &circuit).expect("keygen_pk should not fail");
    end_timer!(timer);

    let timer = start_timer!(|| "create proof");
    let mut transcript = Blake2bWrite::<_, _, Challenge255<_>>::init(vec![]);
    create_proof(
        &params,
        &pk,
        &[circuit.clone()],
        &[&[]],
        OsRng,
        &mut transcript,
    )
    .expect("proof generation should not fail");
    let proof = transcript.finalize();
    end_timer!(timer);

    let timer = start_timer!(|| "gen vk");
    let vk = keygen_vk(&params, &circuit).expect("keygen_vk should not fail");
    end_timer!(timer);

    let timer = start_timer!(|| "verify");
    let params = params.verifier::<Bn256>(0).unwrap();
    let strategy = halo2_proofs::plonk::SingleVerifier::new(&params);
    let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
    halo2_proofs::plonk::verify_proof::<Bn256, _, _, _>(
        &params,
        &vk,
        strategy,
        &[&[]],
        &mut transcript,
    )
    .unwrap();
    end_timer!(timer);
}
