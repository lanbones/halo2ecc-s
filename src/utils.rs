use halo2_proofs::arithmetic::BaseExt;
use num_bigint::BigUint;

pub fn field_to_bn<F: BaseExt>(f: &F) -> BigUint {
    let mut bytes: Vec<u8> = Vec::new();
    f.write(&mut bytes).unwrap();
    BigUint::from_bytes_le(&bytes[..])
}

pub fn bn_to_field<F: BaseExt>(bn: &BigUint) -> F {
    let modulus = field_to_bn(&-F::one()) + 1u64;
    let bn = bn % &modulus;
    let mut bytes = bn.to_bytes_le();
    bytes.resize((modulus.bits() as usize + 7) / 8, 0);
    let mut bytes = &bytes[..];
    F::read(&mut bytes).unwrap()
}
