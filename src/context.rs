use crate::{
    assign::{AssignedValue, Cell, Chip, ValueSchema},
    circuit::{
        base_chip::{BaseChip, FIXED_COLUMNS, MUL_COLUMNS, VAR_COLUMNS},
        range_chip::RangeChip,
    },
    range_info::{RangeClass, RangeInfo, MAX_CHUNKS},
};
use halo2_proofs::{
    arithmetic::{BaseExt, CurveAffine, FieldExt},
    circuit::{AssignedCell, Region},
    plonk::Error,
};
use std::fmt::{Display, Formatter};
use std::marker::PhantomData;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
pub struct Context<W: BaseExt, N: FieldExt> {
    pub records: Arc<Mutex<Records<N>>>,
    pub base_offset: Box<usize>,
    pub range_offset: Box<usize>,
    pub range_info: Option<Arc<RangeInfo<N>>>,
    _mark: PhantomData<W>,
}

impl<W: BaseExt, N: FieldExt> Display for Context<W, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "(range_offset: {}, base_offset: {})",
            self.range_offset, self.base_offset
        )
    }
}

#[repr(transparent)]
pub struct NativeScalarEccContext<C: CurveAffine>(
    pub Context<<C as CurveAffine>::Base, <C as CurveAffine>::ScalarExt>,
);

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    pub fn new() -> Self {
        Self {
            records: Arc::new(Mutex::new(Records::default())),
            base_offset: Box::new(0),
            range_offset: Box::new(0),
            range_info: None,
            _mark: PhantomData,
        }
    }

    pub fn new_with_range_info() -> Self {
        Self {
            records: Arc::new(Mutex::new(Records::default())),
            base_offset: Box::new(0),
            range_offset: Box::new(0),
            range_info: Some(Arc::new(RangeInfo::<N>::new::<W>())),
            _mark: PhantomData,
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct Records<N: FieldExt> {
    pub base_adv_record: Vec<[(Option<N>, bool); VAR_COLUMNS]>,
    pub base_fix_record: Vec<[Option<N>; FIXED_COLUMNS]>,
    pub base_height: usize,

    pub range_adv_record: Vec<(Option<N>, bool)>,
    pub range_fix_record: Vec<[Option<N>; 2]>,
    pub range_height: usize,

    pub permutations: Vec<(Cell, Cell)>,
}

impl<N: FieldExt> Records<N> {
    fn assign_for_max_row(
        &self,
        region: &mut Region<'_, N>,
        base_chip: &BaseChip<N>,
        max_row: usize,
    ) -> Result<bool, Error> {
        let adv_cell = region.assign_advice(
            || "base",
            base_chip.config.base[0],
            max_row,
            || Ok(N::zero()),
        )?;

        let fixed_cell = region.assign_fixed(
            || "base",
            base_chip.config.coeff[0],
            max_row,
            || Ok(N::zero()),
        )?;

        Ok(adv_cell.value().is_none() && fixed_cell.value().is_none())
    }

    fn _assign_to_base_chip(
        &self,
        region: &mut Region<'_, N>,
        base_chip: &BaseChip<N>,
    ) -> Result<Vec<Vec<Option<AssignedCell<N, N>>>>, Error> {
        let mut cells = vec![];

        cells.resize(VAR_COLUMNS, vec![None; self.base_height]);

        for (row, advs) in self.base_adv_record.iter().enumerate() {
            if row >= self.base_height {
                break;
            }

            for (col, adv) in advs.iter().enumerate() {
                if adv.0.is_some() {
                    let cell = region.assign_advice(
                        || "base",
                        base_chip.config.base[col],
                        row,
                        || Ok(adv.0.unwrap()),
                    )?;
                    if adv.1 {
                        cells[col][row] = Some(cell);
                    }
                }
            }
        }

        for (row, fixes) in self.base_fix_record.iter().enumerate() {
            if row >= self.base_height {
                break;
            }

            for (col, fix) in fixes.iter().enumerate() {
                if fix.is_some() {
                    let col = if col < VAR_COLUMNS {
                        base_chip.config.coeff[col]
                    } else if col - VAR_COLUMNS < MUL_COLUMNS {
                        base_chip.config.mul_coeff[col - VAR_COLUMNS]
                    } else if col - VAR_COLUMNS - MUL_COLUMNS == 0 {
                        base_chip.config.next_coeff
                    } else {
                        base_chip.config.constant
                    };

                    region.assign_fixed(|| "fix", col, row, || Ok(fix.unwrap()))?;
                }
            }
        }

        Ok(cells)
    }

    pub fn _assign_to_range_chip<W: BaseExt>(
        &self,
        region: &mut Region<'_, N>,
        range_chip: &RangeChip<W, N>,
    ) -> Result<Vec<Vec<Option<AssignedCell<N, N>>>>, Error> {
        let mut cells = vec![vec![None; self.range_height]];

        for (row, fix) in self.range_fix_record.iter().enumerate() {
            if row >= self.range_height {
                break;
            }
            if fix[0].is_some() {
                region.assign_fixed(
                    || "range block first",
                    range_chip.config.block_first,
                    row,
                    || Ok(fix[0].unwrap()),
                )?;
            }

            if fix[1].is_some() {
                region.assign_fixed(
                    || "range class",
                    range_chip.config.range_class,
                    row,
                    || Ok(fix[1].unwrap()),
                )?;
            }
        }

        for (row, adv) in self.range_adv_record.iter().enumerate() {
            if row >= self.range_height {
                break;
            }
            if adv.0.is_some() {
                let cell = region.assign_advice(
                    || "range var",
                    range_chip.config.value,
                    row,
                    || Ok(adv.0.unwrap()),
                )?;
                if adv.1 {
                    cells[0][row] = Some(cell);
                }
            }
        }

        Ok(cells)
    }

    pub fn _assign_permutation(
        &self,
        region: &mut Region<'_, N>,
        cells: &Vec<Vec<Vec<Option<AssignedCell<N, N>>>>>,
    ) -> Result<(), Error> {
        for (left, right) in self.permutations.iter() {
            let left = cells[left.region as usize][left.col][left.row]
                .as_ref()
                .unwrap()
                .cell();
            let right = cells[right.region as usize][right.col][right.row]
                .as_ref()
                .unwrap()
                .cell();
            region.constrain_equal(left, right)?;
        }

        Ok(())
    }

    pub fn assign_all_in_base(
        &self,
        region: &mut Region<'_, N>,
        base_chip: &BaseChip<N>,
    ) -> Result<Vec<Vec<Vec<Option<AssignedCell<N, N>>>>>, Error> {
        let cells = self._assign_to_base_chip(region, base_chip)?;
        let cells = vec![cells];
        self._assign_permutation(region, &cells)?;
        Ok(cells)
    }

    pub fn assign_all_opt<W: BaseExt>(
        &self,
        region: &mut Region<'_, N>,
        base_chip: &BaseChip<N>,
        range_chip: &RangeChip<W, N>,
    ) -> Result<Option<Vec<Vec<Vec<Option<AssignedCell<N, N>>>>>>, Error> {
        let max_row = usize::max(self.base_height, self.range_height);
        let is_assign_for_max_row = self.assign_for_max_row(region, base_chip, max_row)?;
        if !is_assign_for_max_row {
            let base_cells = self._assign_to_base_chip(region, base_chip)?;
            let range_cells = self._assign_to_range_chip(region, range_chip)?;
            let cells = vec![base_cells, range_cells];
            self._assign_permutation(region, &cells)?;
            Ok(Some(cells))
        } else {
            Ok(None)
        }
    }

    pub fn assign_all<W: BaseExt>(
        &self,
        region: &mut Region<'_, N>,
        base_chip: &BaseChip<N>,
        range_chip: &RangeChip<W, N>,
    ) -> Result<Vec<Vec<Vec<Option<AssignedCell<N, N>>>>>, Error> {
        let base_cells = self._assign_to_base_chip(region, base_chip)?;
        let range_cells = self._assign_to_range_chip(region, range_chip)?;
        let cells = vec![base_cells, range_cells];
        self._assign_permutation(region, &cells)?;
        Ok(cells)
    }

    pub fn enable_permute(&mut self, cell: &Cell) {
        match cell.region {
            Chip::BaseChip => self.base_adv_record[cell.row][cell.col].1 = true,
            Chip::RangeChip => self.range_adv_record[cell.row].1 = true,
        }
    }

    pub fn one_line(
        &mut self,
        offset: usize,
        base_coeff_pairs: Vec<(ValueSchema<N>, N)>,
        constant: Option<N>,
        mul_next_coeffs: (Vec<N>, Option<N>),
    ) {
        assert!(base_coeff_pairs.len() <= VAR_COLUMNS);

        const EXTEND_SIZE: usize = 16;

        if offset >= self.base_adv_record.len() {
            let to_len = (offset + EXTEND_SIZE) & !(EXTEND_SIZE - 1);
            self.base_adv_record
                .resize(to_len, [(None, false); VAR_COLUMNS]);
            self.base_fix_record.resize(to_len, [None; FIXED_COLUMNS]);
        }

        if offset >= self.base_height {
            self.base_height = offset + 1;
        }

        for (i, (base, coeff)) in base_coeff_pairs.into_iter().enumerate() {
            match base.cell() {
                Some(cell) => {
                    let idx = Cell::new(Chip::BaseChip, i, offset);

                    self.base_adv_record[offset][i].1 = true;
                    self.enable_permute(&cell);

                    self.permutations.push((cell, idx));
                }
                _ => {}
            }
            self.base_fix_record[offset][i] = Some(coeff);
            self.base_adv_record[offset][i].0 = Some(base.value());
        }

        let (mul_coeffs, next) = mul_next_coeffs;
        for (i, mul_coeff) in mul_coeffs.into_iter().enumerate() {
            self.base_fix_record[offset][VAR_COLUMNS + i] = Some(mul_coeff);
        }

        if next.is_some() {
            self.base_fix_record[offset][VAR_COLUMNS + MUL_COLUMNS] = next;
        }

        if constant.is_some() {
            self.base_fix_record[offset][VAR_COLUMNS + MUL_COLUMNS + 1] = constant;
        }
    }

    pub fn one_line_with_last(
        &mut self,
        offset: usize,
        base_coeff_pairs: Vec<(ValueSchema<N>, N)>,
        tail: (ValueSchema<N>, N),
        constant: Option<N>,
        mul_next_coeffs: (Vec<N>, Option<N>),
    ) {
        assert!(base_coeff_pairs.len() <= VAR_COLUMNS - 1);

        self.one_line(offset, base_coeff_pairs, constant, mul_next_coeffs);

        let (base, coeff) = tail;

        let i = VAR_COLUMNS - 1;
        match base.cell() {
            Some(cell) => {
                let idx = Cell::new(Chip::BaseChip, i, offset);

                self.base_adv_record[offset][i].1 = true;
                self.enable_permute(&cell);

                self.permutations.push((cell, idx));
            }
            _ => {}
        }
        self.base_fix_record[offset][i] = Some(coeff);
        self.base_adv_record[offset][i].0 = Some(base.value());
    }

    pub fn assign_single_range_value(
        &mut self,
        offset: usize,
        v: N,
        leading_class: RangeClass,
    ) -> AssignedValue<N> {
        const EXTEND_SIZE: usize = 16;

        let end_offset = offset + 1 + MAX_CHUNKS;

        if end_offset >= self.range_adv_record.len() {
            let to_len = (end_offset + EXTEND_SIZE) & !(EXTEND_SIZE - 1);
            self.range_adv_record.resize(to_len, (None, false));
            self.range_fix_record.resize(to_len, [None; 2]);
        }

        if end_offset >= self.range_height {
            self.range_height = end_offset + 1;
        }

        self.range_fix_record[offset][1] = Some(N::from(leading_class as u64));
        self.range_adv_record[offset].0 = Some(v);

        AssignedValue::new(Chip::RangeChip, 0, offset, v)
    }

    pub fn assign_range_value(
        &mut self,
        offset: usize,
        (v, chunks): (N, Vec<N>),
        leading_class: RangeClass,
    ) -> AssignedValue<N> {
        const EXTEND_SIZE: usize = 16;

        let end_offset = offset + 1 + MAX_CHUNKS;

        if end_offset >= self.range_adv_record.len() {
            let to_len = (end_offset + EXTEND_SIZE) & !(EXTEND_SIZE - 1);
            self.range_adv_record.resize(to_len, (None, false));
            self.range_fix_record.resize(to_len, [None; 2]);
        }

        if end_offset >= self.range_height {
            self.range_height = end_offset + 1;
        }

        self.range_fix_record[offset][0] = Some(N::one());
        self.range_adv_record[offset].0 = Some(v);

        // a row placeholder
        self.range_fix_record[offset + MAX_CHUNKS][0] = Some(N::zero());

        for i in 0..chunks.len() - 1 {
            self.range_fix_record[offset + 1 + i][1] = Some(N::from(RangeClass::Common as u64));
        }
        self.range_fix_record[offset + 1 + chunks.len() - 1][1] =
            Some(N::from(leading_class as u64));

        for i in 0..chunks.len() {
            self.range_adv_record[offset + 1 + i].0 = Some(chunks[i]);
        }

        AssignedValue::new(Chip::RangeChip, 0, offset, v)
    }
}
