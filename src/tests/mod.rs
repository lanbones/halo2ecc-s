use halo2_proofs::arithmetic::BaseExt;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

pub mod base_chip;
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
