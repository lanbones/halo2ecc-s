use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::arithmetic::FieldExt;

use super::base_chip::BaseChipOps;
use super::ecc_chip::EccBaseIntegerChipWrapper;
use super::ecc_chip::EccChipScalarOps;
use super::integer_chip::IntegerChipOps;
use crate::assign::AssignedCondition;
use crate::assign::AssignedInteger;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::context::GeneralScalarEccContext;
use crate::pair;
use crate::utils::field_to_bn;

impl<C: CurveAffine, N: FieldExt> EccBaseIntegerChipWrapper<C::Base, N>
    for GeneralScalarEccContext<C, N>
{
    fn base_integer_chip(&mut self) -> &mut dyn IntegerChipOps<C::Base, N> {
        &mut self.base_integer_ctx
    }
}

impl<C: CurveAffine, N: FieldExt> EccChipBaseOps<C, N> for GeneralScalarEccContext<C, N> {}

impl<C: CurveAffine, N: FieldExt> EccChipScalarOps<C, N> for GeneralScalarEccContext<C, N> {
    type AssignedScalar = AssignedInteger<C::Scalar, N>;

    fn decompose_scalar<const WINDOW_SIZE: usize>(
        &mut self,
        s: &Self::AssignedScalar,
    ) -> Vec<[AssignedCondition<N>; WINDOW_SIZE]> {
        let zero = N::zero();
        let one = N::one();
        let two = one + one;
        let two_inv = two.invert().unwrap();

        let s = self.scalar_integer_ctx.reduce(&s);
        let mut bits = vec![];

        for l in s.limbs_le {
            let v = field_to_bn(&l.val);
            let mut rest = l.clone();
            for j in 0..self.scalar_integer_ctx.info.limb_bits {
                let b = self.native_ctx.borrow_mut().assign_bit(v.bit(j).into());
                let v = (rest.val - b.0.val) * two_inv;
                rest = self
                    .native_ctx
                    .borrow_mut()
                    .one_line_with_last(
                        vec![pair!(&rest, -one), pair!(&b.0, one)],
                        pair!(v, two),
                        None,
                        (vec![], None),
                    )
                    .1;
                bits.push(b);
            }

            self.native_ctx.borrow_mut().assert_constant(&rest, zero);
        }

        let padding = bits.len() % WINDOW_SIZE;
        if padding != 0 {
            let zero = self.native_ctx.borrow_mut().assign_constant(zero);
            for _ in 0..padding {
                bits.push(AssignedCondition(zero));
            }
        }

        let mut res = bits.chunks(WINDOW_SIZE)
            .map(|x| Vec::from(x).try_into().unwrap())
            .collect::<Vec<_>>();

        res.reverse();

        res
    }
}
