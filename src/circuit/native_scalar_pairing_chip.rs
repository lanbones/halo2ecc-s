use std::marker::PhantomData;

use halo2_proofs::arithmetic::{BaseExt, CurveAffine, FieldExt};

use crate::{
    assign::{AssignedInteger, AssignedPoint},
    context::{EccContext},
};

pub struct AssignedFq<W: BaseExt, N: FieldExt>(AssignedInteger<W, N>);
pub struct AssignedFq2<W: BaseExt, N: FieldExt>((AssignedFq<W, N>, AssignedFq<W, N>));
pub struct AssignedFq6<W: BaseExt, N: FieldExt>(
    (AssignedFq2<W, N>, AssignedFq2<W, N>, AssignedFq2<W, N>),
);
pub struct AssignedFq12<W: BaseExt, N: FieldExt>((AssignedFq6<W, N>, AssignedFq6<W, N>));

//Todo
pub struct AssignedG1Affine<C: CurveAffine, N: FieldExt>(AssignedPoint<C, N>);
pub struct AssignedG2Affine<C: CurveAffine, N: FieldExt>(
    (AssignedG1Affine<C, N>, AssignedG1Affine<C, N>),
);
pub struct AssignedG2Prepared<C: CurveAffine, N: FieldExt>(PhantomData<(C, N)>);

impl<C: CurveAffine> EccContext<C> {

}
