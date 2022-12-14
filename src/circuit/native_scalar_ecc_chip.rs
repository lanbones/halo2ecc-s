use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::arithmetic::Field;
use halo2_proofs::pairing::group::ff::PrimeField;
use num_integer::Integer;

use super::base_chip::BaseChipOps;
use super::ecc_chip::EccChipScalarOps;
use super::integer_chip::IntegerChipOps;
use crate::assign::AssignedCondition;
use crate::assign::AssignedValue;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::context::NativeScalarEccContext;
use crate::pair;
use crate::utils::{bn_to_field, field_to_bn};

impl<C: CurveAffine> EccChipBaseOps<C, C::Scalar> for NativeScalarEccContext<C> {
    fn base_integer_chip(
        &mut self,
    ) -> &mut dyn IntegerChipOps<<C as CurveAffine>::Base, C::Scalar> {
        &mut self.0
    }
}

impl<C: CurveAffine> EccChipScalarOps<C, C::Scalar> for NativeScalarEccContext<C> {
    type AssignedScalar = AssignedValue<C::ScalarExt>;

    fn decompose_scalar<const WINDOW_SIZE: usize>(
        &mut self,
        s: &AssignedValue<C::ScalarExt>,
    ) -> Vec<[AssignedCondition<C::ScalarExt>; WINDOW_SIZE]> {
        let one = C::ScalarExt::one();
        let two = one + &one;
        let four = two + &two;

        let mut bits = vec![];

        let s_bn = field_to_bn(&s.val);
        let mut v = s.clone();

        for i in 0..<C::ScalarExt as PrimeField>::NUM_BITS as u64 / 2 {
            let b0 = if s_bn.bit(i * 2) {
                C::ScalarExt::one()
            } else {
                C::ScalarExt::zero()
            };
            let b1 = if s_bn.bit(i * 2 + 1) {
                C::ScalarExt::one()
            } else {
                C::ScalarExt::zero()
            };
            let b0 = self.0.ctx.borrow_mut().assign_bit(b0);
            let b1 = self.0.ctx.borrow_mut().assign_bit(b1);
            let v_next: C::ScalarExt = bn_to_field(&(&s_bn >> (i * 2 + 2)));

            let cells = self.0.ctx.borrow_mut().one_line_with_last(
                vec![
                    pair!(v_next.clone(), four),
                    pair!(&b1.0, two),
                    pair!(&b0.0, one),
                ],
                pair!(&v, -one),
                None,
                (vec![], None),
            );

            v = cells.0[0];

            bits.push(b0);
            bits.push(b1);
        }

        if <C::ScalarExt as PrimeField>::NUM_BITS.is_odd() {
            self.0.ctx.borrow_mut().assert_bit(&v);
            bits.push(AssignedCondition(v));
        } else {
            self.0
                .ctx
                .borrow_mut()
                .assert_constant(&v, C::ScalarExt::zero())
        }

        let rem = <C::ScalarExt as PrimeField>::NUM_BITS as usize % WINDOW_SIZE;
        if rem > 0 {
            let zero = self
                .0
                .ctx
                .borrow_mut()
                .assign_constant(C::ScalarExt::zero());
            for _ in 0..WINDOW_SIZE - rem {
                bits.push(AssignedCondition(zero));
            }
        }

        let mut res: Vec<_> = bits
            .chunks(WINDOW_SIZE)
            .into_iter()
            .map(|x| x.try_into().unwrap())
            .collect();
        res.reverse();
        res
    }
}
