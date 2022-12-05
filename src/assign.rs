use std::marker::PhantomData;

use halo2_proofs::arithmetic::{BaseExt, CurveAffine, FieldExt};

use crate::range_info::LIMBS;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum Chip {
    BaseChip = 0,
    RangeChip = 1,
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub struct Cell {
    pub region: Chip,
    pub col: usize,
    pub row: usize,
}

impl Cell {
    pub fn new(region: Chip, col: usize, row: usize) -> Self {
        Self { region, col, row }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct AssignedValue<N: FieldExt> {
    pub cell: Cell,
    pub val: N,
}

#[derive(Debug, Copy, Clone)]
pub struct AssignedInteger<W: BaseExt, N: FieldExt> {
    pub limbs_le: [AssignedValue<N>; LIMBS],
    pub native: AssignedValue<N>,
    pub times: usize,
    _mark: PhantomData<W>,
}

#[derive(Clone, Debug)]
pub struct AssignedCurvature<C: CurveAffine, N: FieldExt>(
    pub AssignedInteger<C::Base, N>,
    pub AssignedCondition<N>,
);

#[derive(Clone, Debug)]
pub struct AssignedPoint<C: CurveAffine, N: FieldExt> {
    pub x: AssignedInteger<C::Base, N>,
    pub y: AssignedInteger<C::Base, N>,
    pub z: AssignedCondition<N>,
}

#[derive(Clone, Debug)]
pub struct AssignedPointWithCurvature<C: CurveAffine, N: FieldExt> {
    pub x: AssignedInteger<C::Base, N>,
    pub y: AssignedInteger<C::Base, N>,
    pub z: AssignedCondition<N>,

    pub curvature: AssignedCurvature<C, N>,
}

impl<C: CurveAffine, N: FieldExt> AssignedPointWithCurvature<C, N> {
    pub fn to_point(&self) -> AssignedPoint<C, N> {
        AssignedPoint::new(self.x, self.y, self.z)
    }
}

impl<W: BaseExt, N: FieldExt> AssignedInteger<W, N> {
    pub fn new(
        limbs_le: [AssignedValue<N>; LIMBS],
        native: AssignedValue<N>,
        times: usize,
    ) -> Self {
        Self {
            limbs_le,
            native,
            times,
            _mark: PhantomData,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct AssignedCondition<N: FieldExt>(pub AssignedValue<N>);

impl<'a, N: FieldExt> AssignedValue<N> {
    pub fn new(region: Chip, col: usize, row: usize, val: N) -> Self {
        Self {
            cell: Cell::new(region, col, row),
            val,
        }
    }
}

impl<C: CurveAffine, N: FieldExt> AssignedPoint<C, N> {
    pub fn new(
        x: AssignedInteger<C::Base, N>,
        y: AssignedInteger<C::Base, N>,
        z: AssignedCondition<N>,
    ) -> Self {
        Self { x, y, z }
    }
}

impl<C: CurveAffine, N: FieldExt> AssignedPointWithCurvature<C, N> {
    pub fn new(
        x: AssignedInteger<C::Base, N>,
        y: AssignedInteger<C::Base, N>,
        z: AssignedCondition<N>,
        curvature: AssignedCurvature<C, N>,
    ) -> Self {
        Self { x, y, z, curvature }
    }
}

#[derive(Debug, Clone)]
pub enum ValueSchema<'a, N: FieldExt> {
    Assigned(&'a AssignedValue<N>),
    Unassigned(N),
}

impl<'a, N: FieldExt> ValueSchema<'a, N> {
    pub fn value(self) -> N {
        match self {
            ValueSchema::Assigned(v) => v.val,
            ValueSchema::Unassigned(v) => v,
        }
    }

    pub fn cell(&self) -> Option<Cell> {
        match self {
            ValueSchema::Assigned(AssignedValue { cell, .. }) => Some(*cell),
            ValueSchema::Unassigned(_) => None,
        }
    }
}

#[macro_export]
macro_rules! pair {
    ($x: expr, $y: expr) => {
        (crate::assign::ValueSchema::from($x), ($y))
    };
}

impl<'a, N: FieldExt> From<N> for ValueSchema<'a, N> {
    fn from(v: N) -> Self {
        Self::Unassigned(v)
    }
}

impl<'a, N: FieldExt> From<&'a AssignedValue<N>> for ValueSchema<'a, N> {
    fn from(v: &'a AssignedValue<N>) -> Self {
        Self::Assigned(v)
    }
}

pub type AssignedFq<W: BaseExt, N: FieldExt> = AssignedInteger<W, N>;
pub type AssignedFq2<W: BaseExt, N: FieldExt> = (AssignedFq<W, N>, AssignedFq<W, N>);
pub type AssignedFq6<W: BaseExt, N: FieldExt> =
    (AssignedFq2<W, N>, AssignedFq2<W, N>, AssignedFq2<W, N>);
pub type AssignedFq12<W: BaseExt, N: FieldExt> = (AssignedFq6<W, N>, AssignedFq6<W, N>);

pub type AssignedG1Affine<C: CurveAffine, N: FieldExt> = AssignedPoint<C, N>;
pub struct AssignedG2Affine<C: CurveAffine, N: FieldExt>(
    pub (AssignedG1Affine<C, N>, AssignedG1Affine<C, N>),
);

//Todo
pub struct AssignedG2Prepared<C: CurveAffine, N: FieldExt>(PhantomData<(C, N)>);
