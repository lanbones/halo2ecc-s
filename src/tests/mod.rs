use halo2_proofs::{
    arithmetic::FieldExt,
    pairing::bn256::{Fq, Fr},
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

pub mod base_chip;
pub mod range_chip;
pub mod integer_chip;
pub mod native_scalar_ecc_chip;
pub mod native_scalar_pairing_chip;

fn random<N: FieldExt>() -> N {
    let seed = chrono::offset::Utc::now()
        .timestamp_nanos()
        .try_into()
        .unwrap();
    let rng = XorShiftRng::seed_from_u64(seed);
    N::random(rng)
}

fn random_fr() -> Fr {
    random::<Fr>()
}

fn random_fq() -> Fq {
    random::<Fq>()
}
