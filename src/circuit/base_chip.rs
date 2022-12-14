use halo2_proofs::{
    arithmetic::FieldExt,
    plonk::{Advice, Column, ConstraintSystem, Fixed},
    poly::Rotation,
};
use std::marker::PhantomData;

use crate::{
    assign::{AssignedCondition, AssignedValue, Chip, ValueSchema},
    context::Context,
    pair,
};

pub const VAR_COLUMNS: usize = 5;
pub const MUL_COLUMNS: usize = 2;
pub const FIXED_COLUMNS: usize = VAR_COLUMNS + MUL_COLUMNS + 2;

#[derive(Clone, Debug)]
pub struct BaseChipConfig {
    pub base: [Column<Advice>; VAR_COLUMNS],
    pub coeff: [Column<Fixed>; VAR_COLUMNS],
    pub mul_coeff: [Column<Fixed>; MUL_COLUMNS],
    pub next_coeff: Column<Fixed>,
    pub constant: Column<Fixed>,
}

#[derive(Clone, Debug)]
pub struct BaseChip<N: FieldExt> {
    pub config: BaseChipConfig,
    mark: PhantomData<N>,
}

impl<N: FieldExt> BaseChip<N> {
    pub fn new(config: BaseChipConfig) -> Self {
        Self {
            config,
            mark: PhantomData,
        }
    }

    pub fn configure(meta: &mut ConstraintSystem<N>) -> BaseChipConfig {
        let base = [(); VAR_COLUMNS].map(|_| meta.advice_column());
        let coeff = [(); VAR_COLUMNS].map(|_| meta.fixed_column());
        let mul_coeff = [(); MUL_COLUMNS].map(|_| meta.fixed_column());
        let next_coeff = meta.fixed_column();
        let constant = meta.fixed_column();

        base.iter().for_each(|c| meta.enable_equality(c.clone()));

        meta.create_gate("base_gate", |meta| {
            let _constant = meta.query_fixed(constant, Rotation::cur());
            let _next = meta.query_advice(base[VAR_COLUMNS - 1], Rotation::next());
            let _next_coeff = meta.query_fixed(next_coeff, Rotation::cur());

            let mut acc = _constant + _next * _next_coeff;
            for i in 0..VAR_COLUMNS {
                let _base = meta.query_advice(base[i], Rotation::cur());
                let _coeff = meta.query_fixed(coeff[i], Rotation::cur());
                acc = acc + _base * _coeff;
            }
            for i in 0..MUL_COLUMNS {
                let _base_l = meta.query_advice(base[i * 2], Rotation::cur());
                let _base_r = meta.query_advice(base[i * 2 + 1], Rotation::cur());
                let _mul_coeff = meta.query_fixed(mul_coeff[i], Rotation::cur());
                acc = acc + _base_l * _base_r * _mul_coeff;
            }

            vec![acc]
        });

        BaseChipConfig {
            base,
            coeff,
            mul_coeff,
            constant,
            next_coeff,
        }
    }
}

pub trait BaseChipOps<N: FieldExt> {
    fn var_columns(&mut self) -> usize;
    fn mul_columns(&mut self) -> usize;

    fn one_line(
        &mut self,
        base_coeff_pairs: Vec<(ValueSchema<'_, N>, N)>,
        constant: Option<N>,
        mul_next_coeffs: (Vec<N>, Option<N>),
    ) -> Vec<AssignedValue<N>>;

    fn one_line_add(
        &mut self,
        base_coeff_pairs: Vec<(ValueSchema<N>, N)>,
        constant: Option<N>,
    ) -> Vec<AssignedValue<N>> {
        self.one_line(base_coeff_pairs, constant, (vec![], None))
    }

    fn one_line_with_last(
        &mut self,
        base_coeff_pairs: Vec<(ValueSchema<'_, N>, N)>,
        last: (ValueSchema<'_, N>, N),
        constant: Option<N>,
        mul_next_coeffs: (Vec<N>, Option<N>),
    ) -> (Vec<AssignedValue<N>>, AssignedValue<N>);

    fn sum_with_constant_in_one_line(
        &mut self,
        elems: Vec<(&AssignedValue<N>, N)>,
        constant: Option<N>,
    ) -> AssignedValue<N> {
        assert!(elems.len() < self.var_columns());

        let sum = elems
            .iter()
            .map(|(x, y)| x.val * y)
            .reduce(|acc, x| acc + x)
            .unwrap();
        let sum = constant.map_or_else(|| sum, |x| x + sum);

        let cells = self.one_line_with_last(
            elems.into_iter().map(|(x, y)| pair!(x, y)).collect(),
            pair!(sum, -N::one()),
            constant,
            (vec![], None),
        );

        cells.1
    }

    fn sum_with_constant(
        &mut self,
        elems: Vec<(&AssignedValue<N>, N)>,
        constant: Option<N>,
    ) -> AssignedValue<N> {
        let columns = self.var_columns();

        if elems.len() < columns {
            self.sum_with_constant_in_one_line(elems, constant)
        } else {
            let (curr, tail) = elems.split_at(columns - 1);
            let mut acc = self.sum_with_constant_in_one_line(Vec::from(curr), constant);

            for chunk in tail.chunks(columns - 2) {
                let elems = vec![chunk, &[(&acc, N::one())]].concat();
                acc = self.sum_with_constant_in_one_line(elems, None);
            }
            acc
        }
    }

    fn add(&mut self, a: &AssignedValue<N>, b: &AssignedValue<N>) -> AssignedValue<N> {
        assert!(self.var_columns() >= 3);

        let one = N::one();
        self.sum_with_constant(vec![(a, one), (b, one)], None)
    }

    fn add_constant(&mut self, a: &AssignedValue<N>, c: N) -> AssignedValue<N> {
        assert!(self.var_columns() >= 3);

        let one = N::one();
        self.sum_with_constant(vec![(a, one)], Some(c))
    }

    fn sub(&mut self, a: &AssignedValue<N>, b: &AssignedValue<N>) -> AssignedValue<N> {
        assert!(self.var_columns() >= 3);

        let one = N::one();
        self.sum_with_constant(vec![(a, one), (b, -one)], None)
    }

    fn mul(&mut self, a: &AssignedValue<N>, b: &AssignedValue<N>) -> AssignedValue<N> {
        assert!(self.var_columns() >= 3);
        assert!(self.mul_columns() >= 1);

        let one = N::one();
        let zero = N::zero();

        let c = a.val * b.val;

        let cells = self.one_line_with_last(
            vec![pair!(a, zero), pair!(b, zero)],
            pair!(c, -one),
            None,
            (vec![one], None),
        );

        cells.1
    }

    fn mul_add_constant(
        &mut self,
        a: &AssignedValue<N>,
        b: &AssignedValue<N>,
        c: N,
    ) -> AssignedValue<N> {
        assert!(self.var_columns() >= 4);
        assert!(self.mul_columns() >= 1);

        let one = N::one();
        let zero = N::zero();

        let d = a.val * b.val + c;

        let cells = self.one_line_with_last(
            vec![pair!(a, zero), pair!(b, zero)],
            pair!(d, -one),
            Some(c),
            (vec![one], None),
        );

        cells.1
    }

    fn mul_add(
        &mut self,
        a: &AssignedValue<N>,
        b: &AssignedValue<N>,
        ab_coeff: N,
        c: &AssignedValue<N>,
        c_coeff: N,
    ) -> AssignedValue<N> {
        assert!(self.var_columns() >= 4);
        assert!(self.mul_columns() >= 1);

        let one = N::one();
        let zero = N::zero();

        let d = a.val * b.val * ab_coeff + c.val * c_coeff;

        let cells = self.one_line_with_last(
            vec![pair!(a, zero), pair!(b, zero), pair!(c, c_coeff)],
            pair!(d, -one),
            None,
            (vec![ab_coeff], None),
        );

        cells.1
    }

    fn mul_add_with_next_line(
        &mut self,
        ls: Vec<(&AssignedValue<N>, &AssignedValue<N>, &AssignedValue<N>, N)>,
    ) -> AssignedValue<N> {
        assert!(self.var_columns() >= 4);
        assert!(self.mul_columns() >= 1);
        assert!(ls.len() > 0);

        if ls.len() == 1 {
            let (a, b, c, c_coeff) = ls[0];
            self.mul_add(a, b, N::one(), c, c_coeff)
        } else {
            let one = N::one();
            let zero = N::zero();

            let mut t = zero;

            for (a, b, c, c_coeff) in ls {
                self.one_line_with_last(
                    vec![pair!(a, zero), pair!(b, zero), pair!(c, c_coeff)],
                    pair!(t, one),
                    None,
                    (vec![one], Some(-one)),
                );

                t = a.val * b.val + c.val * c_coeff + t;
            }

            let cells = self.one_line_with_last(vec![], pair!(t, zero), None, (vec![], None));

            cells.1
        }
    }

    fn invert_unsafe(&mut self, a: &AssignedValue<N>) -> AssignedValue<N> {
        let b = a.val.invert().unwrap();

        let one = N::one();
        let zero = N::zero();

        let cells = self.one_line(
            vec![pair!(a, zero), pair!(b, zero)],
            Some(-one),
            (vec![one], None),
        );

        cells[1]
    }

    fn invert(&mut self, a: &AssignedValue<N>) -> (AssignedCondition<N>, AssignedValue<N>) {
        let zero = N::zero();
        let one = N::one();
        let b = a.val.invert().unwrap_or(zero);
        let c = one - a.val * b;

        // a * c = 0, one of them must be zero
        let cells = self.one_line(
            vec![pair!(a, zero), pair!(c, zero)],
            None,
            (vec![one], None),
        );
        let c = &cells[1];

        // a * b + c = 1
        let cells = self.one_line_with_last(
            vec![pair!(a, zero), pair!(b, zero)],
            pair!(c, one),
            Some(-one),
            (vec![one], None),
        );

        (AssignedCondition(cells.1), cells.0[1])
    }

    fn is_zero(&mut self, a: &AssignedValue<N>) -> AssignedCondition<N> {
        self.invert(a).0
    }

    fn div_unsafe(&mut self, a: &AssignedValue<N>, b: &AssignedValue<N>) -> AssignedValue<N> {
        let c = b.val.invert().unwrap() * a.val;

        let one = N::one();
        let zero = N::zero();

        // b * c = a
        let cells = self.one_line_with_last(
            vec![pair!(b, zero), pair!(c, zero)],
            pair!(a, -one),
            None,
            (vec![one], None),
        );

        cells.0[1]
    }

    fn assign_constant(&mut self, v: N) -> AssignedValue<N> {
        let one = N::one();
        let cells = self.one_line_add(vec![pair!(v, -one)], Some(v));

        cells[0]
    }

    fn assign(&mut self, v: N) -> AssignedValue<N> {
        let zero = N::zero();
        let cells = self.one_line_add(vec![pair!(v, zero)], None);
        cells[0]
    }

    fn assign_bit(&mut self, a: N) -> AssignedCondition<N> {
        let zero = N::zero();
        let one = N::one();

        let cells = self.one_line(
            vec![pair!(a, one), pair!(a, zero)],
            None,
            (vec![-one], None),
        );
        AssignedCondition(cells[0])
    }

    fn assert_equal(&mut self, a: &AssignedValue<N>, b: &AssignedValue<N>) {
        let one = N::one();

        self.one_line_add(vec![pair!(a, -one), pair!(b, one)], None);
    }

    fn assert_constant(&mut self, a: &AssignedValue<N>, b: N) {
        let one = N::one();
        assert_eq!(a.val, b);
        self.one_line_add(vec![pair!(a, -one)], Some(b));
    }

    fn assert_bit(&mut self, a: &AssignedValue<N>) {
        let zero = N::zero();
        let one = N::one();

        self.one_line(
            vec![pair!(a, one), pair!(a, zero)],
            None,
            (vec![-one], None),
        );
    }

    fn and(&mut self, a: &AssignedCondition<N>, b: &AssignedCondition<N>) -> AssignedCondition<N> {
        let res = self.mul(&a.0, &b.0);

        AssignedCondition(res)
    }

    fn not(&mut self, a: &AssignedCondition<N>) -> AssignedCondition<N> {
        let one = N::one();
        let res = self.sum_with_constant(vec![(&a.0, -one)], Some(one));

        AssignedCondition(res)
    }

    fn or(&mut self, a: &AssignedCondition<N>, b: &AssignedCondition<N>) -> AssignedCondition<N> {
        let one = N::one();
        let c = a.0.val + b.0.val - a.0.val * b.0.val;
        let cells = self.one_line_with_last(
            vec![pair!(&a.0, one), pair!(&b.0, one)],
            pair!(c, -one),
            None,
            (vec![-one], None),
        );

        AssignedCondition(cells.1)
    }

    fn xor(&mut self, a: &AssignedCondition<N>, b: &AssignedCondition<N>) -> AssignedCondition<N> {
        let one = N::one();
        let two = one + one;
        let c = a.0.val + b.0.val - two * a.0.val * b.0.val;
        let cells = self.one_line_with_last(
            vec![pair!(&a.0, one), pair!(&b.0, one)],
            pair!(c, -one),
            None,
            (vec![-two], None),
        );

        AssignedCondition(cells.1)
    }

    fn xnor(&mut self, a: &AssignedCondition<N>, b: &AssignedCondition<N>) -> AssignedCondition<N> {
        let one = N::one();
        let two = one + one;
        let c = one - a.0.val - b.0.val + two * a.0.val * b.0.val;
        let cells = self.one_line_with_last(
            vec![pair!(&a.0, -one), pair!(&b.0, -one)],
            pair!(c, -one),
            Some(one),
            (vec![two], None),
        );

        AssignedCondition(cells.1)
    }

    // if cond then a else b
    fn bisec(
        &mut self,
        cond: &AssignedCondition<N>,
        a: &AssignedValue<N>,
        b: &AssignedValue<N>,
    ) -> AssignedValue<N>;

    fn bisec_cond(
        &mut self,
        cond: &AssignedCondition<N>,
        a: &AssignedCondition<N>,
        b: &AssignedCondition<N>,
    ) -> AssignedCondition<N> {
        let c = self.bisec(cond, &a.0, &b.0);
        AssignedCondition(c)
    }

    fn assert_true(&mut self, a: &AssignedCondition<N>) {
        assert!(a.0.val == N::one());
        self.assert_constant(&a.0, N::one())
    }

    fn assert_false(&mut self, a: &AssignedCondition<N>) {
        assert!(a.0.val == N::zero());
        self.assert_constant(&a.0, N::zero())
    }
}

impl<N: FieldExt> BaseChipOps<N> for Context<N> {
    fn var_columns(&mut self) -> usize {
        VAR_COLUMNS
    }

    fn mul_columns(&mut self) -> usize {
        MUL_COLUMNS
    }

    fn one_line(
        &mut self,
        base_coeff_pairs: Vec<(ValueSchema<'_, N>, N)>,
        constant: Option<N>,
        mul_next_coeffs: (Vec<N>, Option<N>),
    ) -> Vec<AssignedValue<N>> {
        let mut records = self.records.lock().unwrap();
        let res = base_coeff_pairs
            .iter()
            .map(|x| x.0.clone().value())
            .enumerate()
            .map(|(i, v)| AssignedValue::new(Chip::BaseChip, i, self.base_offset, v))
            .collect();

        records.one_line(
            self.base_offset,
            base_coeff_pairs,
            constant,
            mul_next_coeffs,
        );

        self.base_offset += 1;

        res
    }

    fn one_line_with_last(
        &mut self,

        base_coeff_pairs: Vec<(ValueSchema<'_, N>, N)>,
        last: (ValueSchema<'_, N>, N),
        constant: Option<N>,
        mul_next_coeffs: (Vec<N>, Option<N>),
    ) -> (Vec<AssignedValue<N>>, AssignedValue<N>) {
        let mut records = self.records.lock().unwrap();
        let res0 = base_coeff_pairs
            .iter()
            .map(|x| x.0.clone().value())
            .enumerate()
            .map(|(i, v)| AssignedValue::new(Chip::BaseChip, i, self.base_offset, v))
            .collect();
        let res1 = AssignedValue::new(
            Chip::BaseChip,
            VAR_COLUMNS - 1,
            self.base_offset,
            last.0.clone().value(),
        );

        records.one_line_with_last(
            self.base_offset,
            base_coeff_pairs,
            last,
            constant,
            mul_next_coeffs,
        );

        self.base_offset += 1;

        (res0, res1)
    }

    fn bisec(
        &mut self,
        cond: &AssignedCondition<N>,
        a: &AssignedValue<N>,
        b: &AssignedValue<N>,
    ) -> AssignedValue<N> {
        let zero = N::zero();
        let one = N::one();

        if self.var_columns() >= 5 {
            let cond_v = cond.0;
            let c = cond.0.val * a.val + (one - cond.0.val) * b.val;
            let cells = self.one_line_with_last(
                vec![
                    pair!(&cond_v, zero),
                    pair!(a, zero),
                    pair!(&cond_v, zero),
                    pair!(b, one),
                ],
                pair!(c, -one),
                None,
                (vec![one, -one], None),
            );

            cells.1
        } else {
            assert!(self.var_columns() >= 3);
            let t = self.mul_add(&cond.0, a, N::one(), b, one);
            self.mul_add(&cond.0, b, -one, &t, one)
        }
    }
}
