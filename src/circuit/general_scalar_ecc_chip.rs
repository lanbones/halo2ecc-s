use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::arithmetic::FieldExt;

use super::ecc_chip::EccBaseIntegerChipWrapper;
use super::ecc_chip::EccChipScalarOps;
use super::integer_chip::IntegerChipOps;
use crate::assign::AssignedCondition;
use crate::assign::AssignedValue;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::context::GeneralScalarEccContext;

impl<C: CurveAffine, N: FieldExt> EccBaseIntegerChipWrapper<C::Base, N>
    for GeneralScalarEccContext<C, N>
{
    fn base_integer_chip(&mut self) -> &mut dyn IntegerChipOps<C::Base, N> {
        &mut self.base_integer_ctx
    }
}

impl<C: CurveAffine, N: FieldExt> EccChipBaseOps<C, N> for GeneralScalarEccContext<C, N> {}

impl<C: CurveAffine, N: FieldExt> EccChipScalarOps<C, N> for GeneralScalarEccContext<C, N> {
    type AssignedScalar = AssignedValue<N>;

    fn decompose_scalar<const WINDOW_SIZE: usize>(
        &mut self,
        _s: &AssignedValue<N>,
    ) -> Vec<[AssignedCondition<N>; WINDOW_SIZE]> {
        todo!()
    }
}
