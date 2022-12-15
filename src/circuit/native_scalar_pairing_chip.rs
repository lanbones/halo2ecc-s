use super::fq12::{Fq12ChipOps, Fq2ChipOps, Fq6ChipOps};
use super::pairing_chip::PairingChipOps;
use crate::context::NativeScalarEccContext;
use halo2_proofs::arithmetic::CurveAffine;

impl<C: CurveAffine> Fq2ChipOps<C::Base, C::Scalar> for NativeScalarEccContext<C> {}
impl<C: CurveAffine> Fq6ChipOps<C::Base, C::Scalar> for NativeScalarEccContext<C> {}
impl<C: CurveAffine> Fq12ChipOps<C::Base, C::Scalar> for NativeScalarEccContext<C> {}
impl<C: CurveAffine> PairingChipOps<C, C::Scalar> for NativeScalarEccContext<C> {}
