use super::fq12::{Fq12ChipOps, Fq2ChipOps, Fq6ChipOps};
use super::pairing_chip::PairingChipOps;
use crate::context::GeneralScalarEccContext;
use halo2_proofs::arithmetic::{CurveAffine, FieldExt};

impl<C: CurveAffine, N: FieldExt> Fq2ChipOps<C::Base, N> for GeneralScalarEccContext<C, N> {}
impl<C: CurveAffine, N: FieldExt> Fq6ChipOps<C::Base, N> for GeneralScalarEccContext<C, N> {}
impl<C: CurveAffine, N: FieldExt> Fq12ChipOps<C::Base, N> for GeneralScalarEccContext<C, N> {}
impl<C: CurveAffine, N: FieldExt> PairingChipOps<C, N> for GeneralScalarEccContext<C, N> {}
