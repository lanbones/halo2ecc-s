/*
  The implementation is ported from https://github.com/privacy-scaling-explorations/pairing
*/
use crate::assign::{AssignedFq6, AssignedG2, AssignedG2Affine, AssignedG2Prepared};
use crate::circuit::integer_chip::IntegerChipOps;
use crate::context::Context;
use crate::utils::bn_to_field;
use crate::{
    assign::{AssignedFq12, AssignedFq2, AssignedG1Affine},
    context::NativeScalarEccContext,
};
use halo2_proofs::arithmetic::{BaseExt, CurveAffine, Engine, FieldExt};
use num_bigint::BigUint;

use super::base_chip::BaseChipOps;

// currently only supports bn256
// 6U+2 for in NAF form
pub const BN_X: u64 = 4965661367192848881;

pub const SIX_U_PLUS_2_NAF: [i8; 65] = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0,
    0, 1, 0, 1, 1,
];

pub const XI_TO_Q_MINUS_1_OVER_2: [[u64; 4]; 2] = [
    [
        0xe4bbdd0c2936b629,
        0xbb30f162e133bacb,
        0x31a9d1b6f9645366,
        0x253570bea500f8dd,
    ],
    [
        0xa1d77ce45ffe77c7,
        0x07affd117826d1db,
        0x6d16bd27bb7edc6b,
        0x2c87200285defecc,
    ],
];

pub const FROBENIUS_COEFF_FQ2_C1: [[u64; 4]; 2] = [
    [
        0xd35d438dc58f0d9d,
        0x0a78eb28f5c70b3d,
        0x666ea36f7879462c,
        0x0e0a77c19a07df2f,
    ],
    [
        0x68c3488912edefaa,
        0x8d087f6872aabf4f,
        0x51e1a24709081231,
        0x2259d6b14729c0fa,
    ],
];

pub const FROBENIUS_COEFF_FQ6_C1: [[[u64; 4]; 2]; 6] = [
    [
        [
            0xd35d438dc58f0d9d,
            0x0a78eb28f5c70b3d,
            0x666ea36f7879462c,
            0x0e0a77c19a07df2f,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0xb5773b104563ab30,
            0x347f91c8a9aa6454,
            0x7a007127242e0991,
            0x1956bcd8118214ec,
        ],
        [
            0x6e849f1ea0aa4757,
            0xaa1c7b6d89f89141,
            0xb6e713cdfae0ca3a,
            0x26694fbb4e82ebc3,
        ],
    ],
    [
        [
            0x3350c88e13e80b9c,
            0x7dce557cdb5e56b9,
            0x6001b4b8b615564a,
            0x2682e617020217e0,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0xc9af22f716ad6bad,
            0xb311782a4aa662b2,
            0x19eeaf64e248c7f4,
            0x20273e77e3439f82,
        ],
        [
            0xacc02860f7ce93ac,
            0x3933d5817ba76b4c,
            0x69e6188b446c8467,
            0x0a46036d4417cc55,
        ],
    ],
    [
        [
            0x71930c11d782e155,
            0xa6bb947cffbe3323,
            0xaa303344d4741444,
            0x2c3b3f0d26594943,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0xf91aba2654e8e3b1,
            0x4771cb2fdc92ce12,
            0xdcb16ae0fc8bdf35,
            0x274aa195cd9d8be4,
        ],
        [
            0x5cfc50ae18811f8b,
            0x4bb28433cb43988c,
            0x4fd35f13c3b56219,
            0x301949bd2fc8883a,
        ],
    ],
];

pub const FROBENIUS_COEFF_FQ6_C2: [[[u64; 4]; 2]; 6] = [
    [
        [
            0xd35d438dc58f0d9d,
            0x0a78eb28f5c70b3d,
            0x666ea36f7879462c,
            0x0e0a77c19a07df2f,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x7361d77f843abe92,
            0xa5bb2bd3273411fb,
            0x9c941f314b3e2399,
            0x15df9cddbb9fd3ec,
        ],
        [
            0x5dddfd154bd8c949,
            0x62cb29a5a4445b60,
            0x37bc870a0c7dd2b9,
            0x24830a9d3171f0fd,
        ],
    ],
    [
        [
            0x71930c11d782e155,
            0xa6bb947cffbe3323,
            0xaa303344d4741444,
            0x2c3b3f0d26594943,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x448a93a57b6762df,
            0xbfd62df528fdeadf,
            0xd858f5d00e9bd47a,
            0x06b03d4d3476ec58,
        ],
        [
            0x2b19daf4bcc936d1,
            0xa1a54e7a56f4299f,
            0xb533eee05adeaef1,
            0x170c812b84dda0b2,
        ],
    ],
    [
        [
            0x3350c88e13e80b9c,
            0x7dce557cdb5e56b9,
            0x6001b4b8b615564a,
            0x2682e617020217e0,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x843420f1d8dadbd6,
            0x31f010c9183fcdb2,
            0x436330b527a76049,
            0x13d47447f11adfe4,
        ],
        [
            0xef494023a857fa74,
            0x2a925d02d5ab101a,
            0x83b015829ba62f10,
            0x2539111d0c13aea3,
        ],
    ],
];

pub const FROBENIUS_COEFF_FQ12_C1: [[[u64; 4]; 2]; 12] = [
    [
        [
            0xd35d438dc58f0d9d,
            0x0a78eb28f5c70b3d,
            0x666ea36f7879462c,
            0x0e0a77c19a07df2f,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0xaf9ba69633144907,
            0xca6b1d7387afb78a,
            0x11bded5ef08a2087,
            0x02f34d751a1f3a7c,
        ],
        [
            0xa222ae234c492d72,
            0xd00f02a4565de15b,
            0xdc2ff3a253dfc926,
            0x10a75716b3899551,
        ],
    ],
    [
        [
            0xca8d800500fa1bf2,
            0xf0c5d61468b39769,
            0x0e201271ad0d4418,
            0x04290f65bad856e6,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x365316184e46d97d,
            0x0af7129ed4c96d9f,
            0x659da72fca1009b5,
            0x08116d8983a20d23,
        ],
        [
            0xb1df4af7c39c1939,
            0x3d9f02878a73bf7f,
            0x9b2220928caf0ae0,
            0x26684515eff054a6,
        ],
    ],
    [
        [
            0x3350c88e13e80b9c,
            0x7dce557cdb5e56b9,
            0x6001b4b8b615564a,
            0x2682e617020217e0,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x86b76f821b329076,
            0x408bf52b4d19b614,
            0x53dfb9d0d985e92d,
            0x051e20146982d2a7,
        ],
        [
            0x0fbc9cd47752ebc7,
            0x6d8fffe33415de24,
            0xbef22cf038cf41b9,
            0x15c0edff3c66bf54,
        ],
    ],
    [
        [
            0x68c3488912edefaa,
            0x8d087f6872aabf4f,
            0x51e1a24709081231,
            0x2259d6b14729c0fa,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x8c84e580a568b440,
            0xcd164d1de0c21302,
            0xa692585790f737d5,
            0x2d7100fdc71265ad,
        ],
        [
            0x99fdddf38c33cfd5,
            0xc77267ed1213e931,
            0xdc2052142da18f36,
            0x1fbcf75c2da80ad7,
        ],
    ],
    [
        [
            0x71930c11d782e155,
            0xa6bb947cffbe3323,
            0xaa303344d4741444,
            0x2c3b3f0d26594943,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0x05cd75fe8a3623ca,
            0x8c8a57f293a85cee,
            0x52b29e86b7714ea8,
            0x2852e0e95d8f9306,
        ],
        [
            0x8a41411f14e0e40e,
            0x59e26809ddfe0b0d,
            0x1d2e2523f4d24d7d,
            0x09fc095cf1414b83,
        ],
    ],
    [
        [
            0x08cfc388c494f1ab,
            0x19b315148d1373d4,
            0x584e90fdcb6c0213,
            0x09e1685bdf2f8849,
        ],
        [0x0, 0x0, 0x0, 0x0],
    ],
    [
        [
            0xb5691c94bd4a6cd1,
            0x56f575661b581478,
            0x64708be5a7fb6f30,
            0x2b462e5e77aecd82,
        ],
        [
            0x2c63ef42612a1180,
            0x29f16aae345bec69,
            0xf95e18c648b216a4,
            0x1aa36073a4cae0d4,
        ],
    ],
];

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    pub fn fq2_assert_equal(&mut self, x: &AssignedFq2<W, N>, y: &AssignedFq2<W, N>) {
        self.assert_int_equal(&x.0, &y.0);
        self.assert_int_equal(&x.1, &y.1);
    }
    pub fn fq2_assign_zero(&mut self) -> AssignedFq2<W, N> {
        let fq2_zero = self.assign_int_constant(W::zero());
        (fq2_zero, fq2_zero)
    }
    pub fn fq2_assign_one(&mut self) -> AssignedFq2<W, N> {
        (
            self.assign_int_constant(W::one()),
            self.assign_int_constant(W::zero()),
        )
    }
    pub fn fq2_assign_constant(&mut self, c0: W, c1: W) -> AssignedFq2<W, N> {
        (self.assign_int_constant(c0), self.assign_int_constant(c1))
    }
    pub fn fq2_add(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_add(&a.0, &b.0), self.int_add(&a.1, &b.1))
    }
    pub fn fq2_mul(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let ab00 = self.int_mul(&a.0, &b.0);
        let ab11 = self.int_mul(&a.1, &b.1);
        let c0 = self.int_sub(&ab00, &ab11);

        let a01 = self.int_add(&a.0, &a.1);
        let b01 = self.int_add(&b.0, &b.1);
        let c1 = self.int_mul(&a01, &b01);
        let c1 = self.int_sub(&c1, &ab00);
        let c1 = self.int_sub(&c1, &ab11);

        (c0, c1)
    }
    pub fn fq2_sub(&mut self, a: &AssignedFq2<W, N>, b: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_sub(&a.0, &b.0), self.int_sub(&a.1, &b.1))
    }
    pub fn fq2_double(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_add(&a.0, &a.0), self.int_add(&a.1, &a.1))
    }
    pub fn fq2_square(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        self.fq2_mul(a, a)
    }
    pub fn fq2_neg(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (self.int_neg(&a.0), self.int_neg(&a.1))
    }
    pub fn fq2_conjugate(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        (a.0, self.int_neg(&a.1))
    }
    pub fn fq2_mul_by_nonresidue(&mut self, a: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let a2 = self.fq2_double(a);
        let a4 = self.fq2_double(&a2);
        let a8 = self.fq2_double(&a4);

        let t = self.int_add(&a8.0, &a.0);
        let c0 = self.int_sub(&t, &a.1);

        let t = self.int_add(&a8.1, &a.0);
        let c1 = self.int_add(&t, &a.1);

        (c0, c1)
    }
    pub fn fq2_unsafe_invert(&mut self, x: &AssignedFq2<W, N>) -> AssignedFq2<W, N> {
        let t0 = self.int_square(&x.0);
        let t1 = self.int_square(&x.1);
        let t0 = self.int_add(&t0, &t1);
        let t = self.int_unsafe_invert(&t0);
        let c0 = self.int_mul(&x.0, &t);
        let c1 = self.int_mul(&x.1, &t);
        let c1 = self.int_neg(&c1);
        (c0, c1)
    }
    pub fn fq2_frobenius_map(&mut self, x: &AssignedFq2<W, N>, power: usize) -> AssignedFq2<W, N> {
        let v = self.assign_int_constant(bn_to_field(&BigUint::from_bytes_le(
            &FROBENIUS_COEFF_FQ2_C1[power % 2]
                .map(|x| x.to_le_bytes())
                .concat(),
        )));
        (x.0, self.int_mul(&x.1, &v))
    }
}

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    pub fn fq6_assert_equal(&mut self, x: &AssignedFq6<W, N>, y: &AssignedFq6<W, N>) {
        self.fq2_assert_equal(&x.0, &y.0);
        self.fq2_assert_equal(&x.1, &y.1);
        self.fq2_assert_equal(&x.2, &y.2);
    }
    pub fn fq6_assign_zero(&mut self) -> AssignedFq6<W, N> {
        let fq2_zero = self.fq2_assign_zero();
        (fq2_zero, fq2_zero, fq2_zero)
    }
    pub fn fq6_assign_one(&mut self) -> AssignedFq6<W, N> {
        let fq2_one = self.fq2_assign_one();
        let fq2_zero = self.fq2_assign_zero();
        (fq2_one, fq2_zero, fq2_zero)
    }
    pub fn fq6_add(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_add(&a.0, &b.0),
            self.fq2_add(&a.1, &b.1),
            self.fq2_add(&a.2, &b.2),
        )
    }
    pub fn fq6_mul(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        let ab00 = self.fq2_mul(&a.0, &b.0);
        let ab11 = self.fq2_mul(&a.1, &b.1);
        let ab22 = self.fq2_mul(&a.2, &b.2);

        let c0 = {
            let b12 = self.fq2_add(&b.1, &b.2);
            let a12 = self.fq2_add(&a.1, &a.2);
            let t = self.fq2_mul(&a12, &b12);
            let t = self.fq2_sub(&t, &ab11);
            let t = self.fq2_sub(&t, &ab22);
            let t = self.fq2_mul_by_nonresidue(&t);
            self.fq2_add(&t, &ab00)
        };

        let c1 = {
            let b01 = self.fq2_add(&b.0, &b.1);
            let a01 = self.fq2_add(&a.0, &a.1);
            let t = self.fq2_mul(&a01, &b01);
            let t = self.fq2_sub(&t, &ab00);
            let t = self.fq2_sub(&t, &ab11);
            let t = self.fq2_mul_by_nonresidue(&t);
            self.fq2_add(&t, &ab22)
        };

        let c2 = {
            let b02 = self.fq2_add(&b.0, &b.2);
            let a02 = self.fq2_add(&a.0, &a.2);
            let t = self.fq2_mul(&a02, &b02);
            let t = self.fq2_sub(&t, &ab00);
            let t = self.fq2_add(&t, &ab11);
            self.fq2_sub(&t, &ab22)
        };

        (c0, c1, c2)
    }
    pub fn fq6_sub(&mut self, a: &AssignedFq6<W, N>, b: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_sub(&a.0, &b.0),
            self.fq2_sub(&a.1, &b.1),
            self.fq2_sub(&a.2, &b.2),
        )
    }
    pub fn fq6_double(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (
            self.fq2_double(&a.0),
            self.fq2_double(&a.1),
            self.fq2_double(&a.2),
        )
    }
    pub fn fq6_square(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        self.fq6_mul(a, a)
    }
    pub fn fq6_neg(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (self.fq2_neg(&a.0), self.fq2_neg(&a.1), self.fq2_neg(&a.2))
    }
    pub fn fq6_mul_by_nonresidue(&mut self, a: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        (self.fq2_mul_by_nonresidue(&a.2), a.0, a.1)
    }
    pub fn fq6_mul_by_1(
        &mut self,
        a: &AssignedFq6<W, N>,
        b1: &AssignedFq2<W, N>,
    ) -> AssignedFq6<W, N> {
        let ab11 = self.fq2_mul(&a.1, &b1);

        let c0 = {
            let b12 = b1;
            let a12 = self.fq2_add(&a.1, &a.2);
            let t = self.fq2_mul(&a12, &b12);
            let t = self.fq2_sub(&t, &ab11);
            self.fq2_mul_by_nonresidue(&t)
        };

        let c1 = {
            let b01 = b1;
            let a01 = self.fq2_add(&a.0, &a.1);
            let t = self.fq2_mul(&a01, &b01);
            self.fq2_sub(&t, &ab11)
        };

        let c2 = ab11;

        (c0, c1, c2)
    }
    pub fn fq6_mul_by_01(
        &mut self,
        a: &AssignedFq6<W, N>,
        b0: &AssignedFq2<W, N>,
        b1: &AssignedFq2<W, N>,
    ) -> AssignedFq6<W, N> {
        let ab00 = self.fq2_mul(&a.0, &b0);
        let ab11 = self.fq2_mul(&a.1, &b1);

        let c0 = {
            let b12 = b1;
            let a12 = self.fq2_add(&a.1, &a.2);
            let t = self.fq2_mul(&a12, &b12);
            let t = self.fq2_sub(&t, &ab11);
            let t = self.fq2_mul_by_nonresidue(&t);
            self.fq2_add(&t, &ab00)
        };

        let c1 = {
            let b01 = self.fq2_add(b0, b1);
            let a01 = self.fq2_add(&a.0, &a.1);
            let t = self.fq2_mul(&a01, &b01);
            let t = self.fq2_sub(&t, &ab00);
            self.fq2_sub(&t, &ab11)
        };

        let c2 = {
            let b02 = b0;
            let a02 = self.fq2_add(&a.0, &a.2);
            let t = self.fq2_mul(&a02, &b02);
            let t = self.fq2_sub(&t, &ab00);
            self.fq2_add(&t, &ab11)
        };

        (c0, c1, c2)
    }
    pub fn fq6_unsafe_invert(&mut self, x: &AssignedFq6<W, N>) -> AssignedFq6<W, N> {
        let c0 = self.fq2_mul_by_nonresidue(&x.2);
        let c0 = self.fq2_mul(&c0, &x.1);
        let c0 = self.fq2_neg(&c0);
        let x0s = self.fq2_square(&x.0);
        let c0 = self.fq2_add(&c0, &x0s);

        let c1 = self.fq2_square(&x.2);
        let c1 = self.fq2_mul_by_nonresidue(&c1);
        let x01 = self.fq2_mul(&x.0, &x.1);
        let c1 = self.fq2_sub(&c1, &x01);

        let c2 = self.fq2_square(&x.1);
        let x02 = self.fq2_mul(&x.0, &x.2);
        let c2 = self.fq2_sub(&c2, &x02);

        let c0x0 = self.fq2_mul(&c0, &x.0);
        let c1x2 = self.fq2_mul(&c1, &x.2);
        let c2x1 = self.fq2_mul(&c2, &x.1);
        let t = self.fq2_add(&c1x2, &c2x1);
        let t = self.fq2_mul_by_nonresidue(&t);
        let t = self.fq2_add(&t, &c0x0);

        let t = self.fq2_unsafe_invert(&t);

        (
            self.fq2_mul(&t, &c0),
            self.fq2_mul(&t, &c1),
            self.fq2_mul(&t, &c2),
        )
    }
    pub fn fq6_frobenius_map(&mut self, x: &AssignedFq6<W, N>, power: usize) -> AssignedFq6<W, N> {
        let c0 = self.fq2_frobenius_map(&x.0, power);
        let c1 = self.fq2_frobenius_map(&x.1, power);
        let c2 = self.fq2_frobenius_map(&x.2, power);

        let coeff_c1 = FROBENIUS_COEFF_FQ6_C1[power % 6].map(|x| {
            bn_to_field(&BigUint::from_bytes_le(
                &x.map(|x| x.to_le_bytes()).concat(),
            ))
        });
        let coeff_c1 = self.fq2_assign_constant(coeff_c1[0], coeff_c1[1]);
        let c1 = self.fq2_mul(&c1, &coeff_c1);
        let coeff_c2 = FROBENIUS_COEFF_FQ6_C2[power % 6].map(|x| {
            bn_to_field(&BigUint::from_bytes_le(
                &x.map(|x| x.to_le_bytes()).concat(),
            ))
        });
        let coeff_c2 = self.fq2_assign_constant(coeff_c2[0], coeff_c2[1]);
        let c2 = self.fq2_mul(&c2, &coeff_c2);

        (c0, c1, c2)
    }
}

impl<W: BaseExt, N: FieldExt> Context<W, N> {
    pub fn fq12_assert_identity(&mut self, x: &AssignedFq12<W, N>) {
        let one = self.fq12_assign_one();
        self.fq12_assert_eq(x, &one);
    }
    pub fn fq12_assert_eq(&mut self, x: &AssignedFq12<W, N>, y: &AssignedFq12<W, N>) {
        self.fq6_assert_equal(&x.0, &y.0);
        self.fq6_assert_equal(&x.1, &y.1);
    }
    pub fn fq12_assign_zero(&mut self) -> AssignedFq12<W, N> {
        let fq6_zero = self.fq6_assign_zero();
        (fq6_zero, fq6_zero)
    }
    pub fn fq12_assign_one(&mut self) -> AssignedFq12<W, N> {
        let fq6_one = self.fq6_assign_one();
        let fq6_zero = self.fq6_assign_zero();
        (fq6_one, fq6_zero)
    }
    pub fn fq12_add(
        &mut self,
        a: &AssignedFq12<W, N>,
        b: &AssignedFq12<W, N>,
    ) -> AssignedFq12<W, N> {
        (self.fq6_add(&a.0, &b.0), self.fq6_add(&a.1, &b.1))
    }
    pub fn fq12_mul(
        &mut self,
        a: &AssignedFq12<W, N>,
        b: &AssignedFq12<W, N>,
    ) -> AssignedFq12<W, N> {
        let ab00 = self.fq6_mul(&a.0, &b.0);
        let ab11 = self.fq6_mul(&a.1, &b.1);

        let a01 = self.fq6_add(&a.0, &a.1);
        let b01 = self.fq6_add(&b.0, &b.1);
        let c1 = self.fq6_mul(&a01, &b01);
        let c1 = self.fq6_sub(&c1, &ab00);
        let c1 = self.fq6_sub(&c1, &ab11);

        let ab11 = self.fq6_mul_by_nonresidue(&ab11);
        let c0 = self.fq6_add(&ab00, &ab11);

        (c0, c1)
    }
    pub fn fq12_sub(
        &mut self,
        a: &AssignedFq12<W, N>,
        b: &AssignedFq12<W, N>,
    ) -> AssignedFq12<W, N> {
        (self.fq6_sub(&a.0, &b.0), self.fq6_sub(&a.1, &b.1))
    }
    pub fn fq12_double(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_double(&a.0), self.fq6_double(&a.1))
    }
    pub fn fq12_square(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        self.fq12_mul(a, a)
    }
    pub fn fq12_neg(&mut self, a: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (self.fq6_neg(&a.0), self.fq6_neg(&a.1))
    }
    pub fn fq12_conjugate(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        (x.0, self.fq6_neg(&x.1))
    }
    pub fn fq12_mul_by_014(
        &mut self,
        x: &AssignedFq12<W, N>,
        c0: &AssignedFq2<W, N>,
        c1: &AssignedFq2<W, N>,
        c4: &AssignedFq2<W, N>,
    ) -> AssignedFq12<W, N> {
        let t0 = self.fq6_mul_by_01(&x.0, c0, c1);
        let t1 = self.fq6_mul_by_1(&x.1, c4);
        let o = self.fq2_add(c1, c4);

        let x0 = self.fq6_add(&x.0, &x.1);
        let x0 = self.fq6_mul_by_01(&x0, c0, &o);
        let x0 = self.fq6_sub(&x0, &t0);
        let x0 = self.fq6_sub(&x0, &t1);

        let x1 = self.fq6_mul_by_nonresidue(&t1);
        let x1 = self.fq6_add(&x1, &t0);

        (x0, x1)
    }
    pub fn fq12_mul_by_034(
        &mut self,
        x: &AssignedFq12<W, N>,
        c0: &AssignedFq2<W, N>,
        c3: &AssignedFq2<W, N>,
        c4: &AssignedFq2<W, N>,
    ) -> AssignedFq12<W, N> {
        let t00 = self.fq2_mul(&x.0 .0, c0);
        let t01 = self.fq2_mul(&x.0 .1, c0);
        let t02 = self.fq2_mul(&x.0 .2, c0);
        let t0 = (t00, t01, t02);

        let t1 = self.fq6_mul_by_01(&x.1, c3, c4);
        let t2 = self.fq6_add(&x.0, &x.1);
        let o = self.fq2_add(c0, c3);
        let t2 = self.fq6_mul_by_01(&t2, &o, c4);
        let t2 = self.fq6_sub(&t2, &t0);
        let x0 = self.fq6_sub(&t2, &t1);
        let t1 = self.fq6_mul_by_nonresidue(&t1);
        let x1 = self.fq6_add(&t0, &t1);
        (x0, x1)
    }
    pub fn fp4_square(
        &mut self,
        c0: &mut AssignedFq2<W, N>,
        c1: &mut AssignedFq2<W, N>,
        a0: &AssignedFq2<W, N>,
        a1: &AssignedFq2<W, N>,
    ) {
        let t0 = self.fq2_square(&a0);
        let t1 = self.fq2_square(&a1);
        let mut t2 = self.fq2_mul_by_nonresidue(&t1);
        *c0 = self.fq2_add(&t2, &t0);
        t2 = self.fq2_add(a0, a1);
        t2 = self.fq2_square(&t2);
        t2 = self.fq2_sub(&t2, &t0);
        *c1 = self.fq2_sub(&t2, &t1);
    }
    pub fn fq12_cyclotomic_square(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        let zero = self.fq2_assign_zero();
        let mut t3 = zero;
        let mut t4 = zero;
        let mut t5 = zero;
        let mut t6 = zero;

        self.fp4_square(&mut t3, &mut t4, &x.0 .0, &x.1 .1);
        let mut t2 = self.fq2_sub(&t3, &x.0 .0);
        t2 = self.fq2_double(&t2);
        let c00 = self.fq2_add(&t2, &t3);

        t2 = self.fq2_add(&t4, &x.1 .1);
        t2 = self.fq2_double(&t2);
        let c11 = self.fq2_add(&t2, &t4);

        self.fp4_square(&mut t3, &mut t4, &x.1 .0, &x.0 .2);
        self.fp4_square(&mut t5, &mut t6, &x.0 .1, &x.1 .2);

        t2 = self.fq2_sub(&t3, &x.0 .1);
        t2 = self.fq2_double(&t2);
        let c01 = self.fq2_add(&t2, &t3);
        t2 = self.fq2_add(&t4, &x.1 .2);
        t2 = self.fq2_double(&t2);
        let c12 = self.fq2_add(&t2, &t4);
        t3 = t6;
        t3 = self.fq2_mul_by_nonresidue(&t3);
        t2 = self.fq2_add(&t3, &x.1 .0);
        t2 = self.fq2_double(&t2);
        let c10 = self.fq2_add(&t2, &t3);
        t2 = self.fq2_sub(&t5, &x.0 .2);
        t2 = self.fq2_double(&t2);
        let c02 = self.fq2_add(&t2, &t5);

        ((c00, c01, c02), (c10, c11, c12))
    }
    pub fn fq12_unsafe_invert(&mut self, x: &AssignedFq12<W, N>) -> AssignedFq12<W, N> {
        let x0s = self.fq6_square(&x.0);
        let x1s = self.fq6_square(&x.1);
        let t = self.fq6_mul_by_nonresidue(&x1s);
        let t = self.fq6_sub(&x0s, &t);
        let t = self.fq6_unsafe_invert(&t);

        let c0 = self.fq6_mul(&t, &x.0);
        let c1 = self.fq6_mul(&t, &x.1);
        let c1 = self.fq6_neg(&c1);
        (c0, c1)
    }
    pub fn fq12_frobenius_map(
        &mut self,
        x: &AssignedFq12<W, N>,
        power: usize,
    ) -> AssignedFq12<W, N> {
        let c0 = self.fq6_frobenius_map(&x.0, power);
        let c1 = self.fq6_frobenius_map(&x.1, power);

        let coeff = FROBENIUS_COEFF_FQ12_C1[power % 12].map(|x| {
            bn_to_field(&BigUint::from_bytes_le(
                &x.map(|x| x.to_le_bytes()).concat(),
            ))
        });
        let coeff = self.fq2_assign_constant(coeff[0], coeff[1]);
        let c1c0 = self.fq2_mul(&c1.0, &coeff);
        let c1c1 = self.fq2_mul(&c1.1, &coeff);
        let c1c2 = self.fq2_mul(&c1.2, &coeff);

        (c0, (c1c0, c1c1, c1c2))
    }
}

impl<C: CurveAffine> NativeScalarEccContext<C> {
    pub fn doubling_step(
        &mut self,
        pt: &AssignedG2<C, C::Scalar>,
    ) -> (
        AssignedG2<C, C::Scalar>,
        [AssignedFq2<C::Base, C::Scalar>; 3],
    ) {
        //see https://eprint.iacr.org/2010/354.pdf
        let x2 = self.0.fq2_square(&pt.x);

        let y2 = self.0.fq2_square(&pt.y);
        let _2y2 = self.0.fq2_double(&y2);
        let _4y2 = self.0.fq2_double(&_2y2);
        let _4y4 = self.0.fq2_square(&_2y2);
        let _8y4 = self.0.fq2_double(&_4y4);

        let z2 = self.0.fq2_square(&pt.z);

        let _4xy2 = {
            let t = self.0.fq2_mul(&y2, &pt.x);
            let t = self.0.fq2_double(&t);
            let t = self.0.fq2_double(&t);
            t
        };
        let _3x2 = {
            let t = self.0.fq2_double(&x2);
            let t = self.0.fq2_add(&t, &x2);
            t
        };
        let _6x2 = self.0.fq2_double(&_3x2);
        let _9x4 = self.0.fq2_square(&_3x2);
        let _3x2_x = self.0.fq2_add(&_3x2, &pt.x);

        let rx = {
            let t = self.0.fq2_sub(&_9x4, &_4xy2);
            let t = self.0.fq2_sub(&t, &_4xy2);
            t
        };

        let ry = {
            let t = self.0.fq2_sub(&_4xy2, &pt.x);
            let t = self.0.fq2_mul(&t, &_3x2_x);
            let t = self.0.fq2_sub(&t, &_8y4);
            t
        };

        let rz = {
            let yz = self.0.fq2_mul(&pt.y, &pt.z);
            self.0.fq2_double(&yz)
        };

        let c0 = {
            let t = self.0.fq2_mul(&z2, &rz);
            self.0.fq2_double(&t)
        };

        let c1 = {
            let _6x2z2 = self.0.fq2_mul(&z2, &_6x2);
            self.0.fq2_neg(&_6x2z2)
        };

        let c2 = {
            let _6x3 = self.0.fq2_mul(&_6x2, &pt.x);
            self.0.fq2_sub(&_6x3, &_4y2)
        };

        (AssignedG2::new(rx, ry, rz), [c0, c1, c2])
    }

    pub fn addition_step(
        &mut self,
        pt: &AssignedG2<C, C::Scalar>,
        pq: &AssignedG2Affine<C, C::Scalar>,
    ) -> (
        AssignedG2<C, C::Scalar>,
        [AssignedFq2<C::Base, C::Scalar>; 3],
    ) {
        //see https://eprint.iacr.org/2010/354.pdf
        let zt2 = self.0.fq2_square(&pt.z);
        let yqzt = self.0.fq2_mul(&pq.y, &pt.z);
        let yqzt3 = self.0.fq2_mul(&yqzt, &zt2);
        let _2yqzt3 = self.0.fq2_double(&yqzt3);
        let _2yt = self.0.fq2_double(&pt.y);
        let _2yqzt3_2yt = self.0.fq2_sub(&_2yqzt3, &_2yt);

        let xqzt2 = self.0.fq2_mul(&pq.x, &zt2);
        let xqzt2_xt = self.0.fq2_sub(&xqzt2, &pt.x);
        let _2_xqzt2_xt = self.0.fq2_double(&xqzt2_xt); // 2(xqzt2 - xt)
        let _4_xqzt2_xt_2 = self.0.fq2_square(&_2_xqzt2_xt); // 4(xqzt2 - xt) ^ 2

        let rx = {
            let t0 = self.0.fq2_mul(&_4_xqzt2_xt_2, &xqzt2_xt); // 4(xqzt2 - xt) ^ 3
            let t1 = self.0.fq2_double(&_4_xqzt2_xt_2); // 8(xqzt2 - xt) ^ 2
            let t2 = self.0.fq2_mul(&t1, &pt.x); // 8(xqzt2 - xt) ^ 2 * xt

            let t = self.0.fq2_square(&_2yqzt3_2yt); // (2yqzt3 - 2yt) ^ 2
            let t = self.0.fq2_sub(&t, &t0); // (2yqzt3 - 2yt) ^ 2 - 4(xqzt2 - xt) ^ 3
            let t = self.0.fq2_sub(&t, &t2); // (2yqzt3 - 2yt) ^ 2 - 4(xqzt2 - xt) ^ 3 - 8(xqzt2 - xt) ^ 2 * xt
            t
        };

        let ry = {
            let t0 = self.0.fq2_mul(&_4_xqzt2_xt_2, &pt.x); // 4(xqzt2 - xt) ^ 2 * xt
            let t0 = self.0.fq2_sub(&t0, &rx); // 4(xqzt2 - xt) ^ 2 * xt - xr
            let t0 = self.0.fq2_mul(&_2yqzt3_2yt, &t0); // (2yqzt3 - 2yt) * (4(xqzt2 - xt) ^ 2 * xt - xr)
            let t1 = self.0.fq2_mul(&_2_xqzt2_xt, &_4_xqzt2_xt_2); // 8(xqzt2 - xt) ^ 3
            let t1 = self.0.fq2_mul(&t1, &pt.y); // 8yt * (xqzt2 - xt) ^ 3
            let t = self.0.fq2_sub(&t0, &t1);
            t
        };

        let rz = self.0.fq2_mul(&pt.z, &_2_xqzt2_xt);

        let c0 = self.0.fq2_double(&rz);
        let c1 = {
            let t = self.0.fq2_add(&_2yqzt3, &_2yt); // 2(yqzt3 + yt)
            self.0.fq2_double(&t)
        };
        let c2 = {
            let t0 = self.0.fq2_mul(&_2yqzt3, &pq.x); // 2yqzt3xq
            let t0 = self.0.fq2_sub(&t0, &_2yt); // 2(yqzt3xq - yt)
            let t0 = self.0.fq2_double(&t0); // 4(yqzt3xq - yt)
            let t0 = self.0.fq2_mul(&t0, &pq.x); // 4xq(yqzt3xq - yt)
            let t1 = self.0.fq2_mul(&pq.y, &rz); // yqzr
            let t1 = self.0.fq2_double(&t1); // 2yqzr
            self.0.fq2_sub(&t0, &t1) // 4xq(yqzt3xq - yt) - yqzr
        };

        (AssignedG2::new(rx, ry, rz), [c0, c1, c2])
    }

    pub fn g2affine_to_g2(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2<C, C::Scalar> {
        // not support identity
        self.0.assert_false(&g2.z);
        let z = self.0.fq2_assign_one();

        AssignedG2::new(g2.x, g2.y, z)
    }

    pub fn g2_neg(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2Affine<C, C::Scalar> {
        let y = self.0.fq2_neg(&g2.y);
        AssignedG2Affine::new(g2.x, y, g2.z)
    }

    pub fn prepare_g2(
        &mut self,
        g2: &AssignedG2Affine<C, C::Scalar>,
    ) -> AssignedG2Prepared<C, C::ScalarExt> {
        let neg_g2 = self.g2_neg(&g2);

        let mut coeffs = vec![];
        let mut r = self.g2affine_to_g2(g2);

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            let (new_r, coeff) = self.doubling_step(&r);
            r = new_r;
            coeffs.push(coeff);
            let x = SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    let (new_r, coeff) = self.addition_step(&r, &g2);
                    r = new_r;
                    coeffs.push(coeff);
                }
                -1 => {
                    let (new_r, coeff) = self.addition_step(&r, &neg_g2);
                    r = new_r;
                    coeffs.push(coeff);
                }
                _ => continue,
            }
        }

        let mut q1 = g2.clone();

        let c11 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(
                &FROBENIUS_COEFF_FQ6_C1[1][0]
                    .map(|x| x.to_le_bytes())
                    .concat(),
            )),
            bn_to_field(&BigUint::from_bytes_le(
                &FROBENIUS_COEFF_FQ6_C1[1][1]
                    .map(|x| x.to_le_bytes())
                    .concat(),
            )),
        );
        let c12 = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(
                &FROBENIUS_COEFF_FQ6_C1[2][0]
                    .map(|x| x.to_le_bytes())
                    .concat(),
            )),
            bn_to_field(&BigUint::from_bytes_le(
                &FROBENIUS_COEFF_FQ6_C1[2][1]
                    .map(|x| x.to_le_bytes())
                    .concat(),
            )),
        );
        let xi = self.0.fq2_assign_constant(
            bn_to_field(&BigUint::from_bytes_le(
                &XI_TO_Q_MINUS_1_OVER_2[0].map(|x| x.to_le_bytes()).concat(),
            )),
            bn_to_field(&BigUint::from_bytes_le(
                &XI_TO_Q_MINUS_1_OVER_2[1].map(|x| x.to_le_bytes()).concat(),
            )),
        );

        q1.x.1 = self.0.int_neg(&q1.x.1);
        q1.x = self.0.fq2_mul(&q1.x, &c11);

        q1.y.1 = self.0.int_neg(&q1.y.1);
        q1.y = self.0.fq2_mul(&q1.y, &xi);

        let (new_r, coeff) = self.addition_step(&r, &q1);
        r = new_r;
        coeffs.push(coeff);

        let mut minusq2 = g2.clone();
        minusq2.x = self.0.fq2_mul(&minusq2.x, &c12);

        let (new_r, coeff) = self.addition_step(&r, &minusq2);
        r = new_r;
        coeffs.push(coeff);

        AssignedG2Prepared::new(coeffs)
    }

    pub fn ell(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
        coeffs: &[AssignedFq2<C::Base, C::Scalar>; 3],
        p: &AssignedG1Affine<C, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let c00 = &coeffs[0].0;
        let c01 = &coeffs[0].1;
        let c10 = &coeffs[1].0;
        let c11 = &coeffs[1].1;

        let c00 = self.0.int_mul(&c00, &p.y);
        let c01 = self.0.int_mul(&c01, &p.y);
        let c10 = self.0.int_mul(&c10, &p.x);
        let c11 = self.0.int_mul(&c11, &p.x);

        self.0
            .fq12_mul_by_034(f, &(c00, c01), &(c10, c11), &coeffs[2])
    }

    pub fn multi_miller_loop(
        &mut self,
        terms: &[(
            &AssignedG1Affine<C, C::Scalar>,
            &AssignedG2Prepared<C, C::Scalar>,
        )],
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let mut pairs = vec![];
        for &(p, q) in terms {
            // not support identity
            self.0.assert_false(&p.z);
            pairs.push((p, q.coeffs.iter()));
        }

        let mut f = self.0.fq12_assign_one();

        for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
            if i != SIX_U_PLUS_2_NAF.len() - 1 {
                f = self.0.fq12_square(&f);
            }
            for &mut (p, ref mut coeffs) in &mut pairs {
                f = self.ell(&f, coeffs.next().unwrap(), &p);
            }
            let x = SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    for &mut (p, ref mut coeffs) in &mut pairs {
                        f = self.ell(&f, coeffs.next().unwrap(), &p);
                    }
                }
                -1 => {
                    for &mut (p, ref mut coeffs) in &mut pairs {
                        f = self.ell(&f, coeffs.next().unwrap(), &p);
                    }
                }
                _ => continue,
            }
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            f = self.ell(&f, coeffs.next().unwrap(), &p);
        }

        for &mut (p, ref mut coeffs) in &mut pairs {
            f = self.ell(&f, coeffs.next().unwrap(), &p);
        }

        for &mut (_p, ref mut coeffs) in &mut pairs {
            assert!(coeffs.next().is_none());
        }

        f
    }

    pub fn exp_by_x(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let x = BN_X;
        let mut res = self.0.fq12_assign_one();
        for i in (0..64).rev() {
            res = self.0.fq12_cyclotomic_square(&res);
            if ((x >> i) & 1) == 1 {
                res = self.0.fq12_mul(&res, &f);
            }
        }
        res
    }

    pub fn final_exponentiation(
        &mut self,
        f: &AssignedFq12<C::Base, C::Scalar>,
    ) -> AssignedFq12<C::Base, C::Scalar> {
        let f1 = self.0.fq12_conjugate(&f);
        let mut f2 = self.0.fq12_unsafe_invert(&f);

        let mut r = self.0.fq12_mul(&f1, &f2);
        f2 = r;
        r = self.0.fq12_frobenius_map(&r, 2);
        r = self.0.fq12_mul(&r, &f2);

        let mut fp = r;
        fp = self.0.fq12_frobenius_map(&fp, 1);

        let mut fp2 = r;
        fp2 = self.0.fq12_frobenius_map(&fp2, 2);
        let mut fp3 = fp2;
        fp3 = self.0.fq12_frobenius_map(&fp3, 1);

        let mut fu = r;
        fu = self.exp_by_x(&fu);

        let mut fu2 = fu;
        fu2 = self.exp_by_x(&fu2);

        let mut fu3 = fu2;
        fu3 = self.exp_by_x(&fu3);

        let mut y3 = fu;
        y3 = self.0.fq12_frobenius_map(&y3, 1);

        let mut fu2p = fu2;
        fu2p = self.0.fq12_frobenius_map(&fu2p, 1);

        let mut fu3p = fu3;
        fu3p = self.0.fq12_frobenius_map(&fu3p, 1);

        let mut y2 = fu2;
        y2 = self.0.fq12_frobenius_map(&y2, 2);

        let mut y0 = fp;
        y0 = self.0.fq12_mul(&y0, &fp2);
        y0 = self.0.fq12_mul(&y0, &fp3);

        let mut y1 = r;
        y1 = self.0.fq12_conjugate(&y1);

        let mut y5 = fu2;
        y5 = self.0.fq12_conjugate(&y5);

        y3 = self.0.fq12_conjugate(&y3);

        let mut y4 = fu;
        y4 = self.0.fq12_mul(&y4, &fu2p);
        y4 = self.0.fq12_conjugate(&y4);

        let mut y6 = fu3;
        y6 = self.0.fq12_mul(&y6, &fu3p);
        y6 = self.0.fq12_conjugate(&y6);

        y6 = self.0.fq12_cyclotomic_square(&y6);
        y6 = self.0.fq12_mul(&y6, &y4);
        y6 = self.0.fq12_mul(&y6, &y5);

        let mut t1 = y3;
        t1 = self.0.fq12_mul(&t1, &y5);
        t1 = self.0.fq12_mul(&t1, &y6);

        y6 = self.0.fq12_mul(&y6, &y2);

        t1 = self.0.fq12_cyclotomic_square(&t1);
        t1 = self.0.fq12_mul(&t1, &y6);
        t1 = self.0.fq12_cyclotomic_square(&t1);

        let mut t0 = t1;
        t0 = self.0.fq12_mul(&t0, &y1);

        t1 = self.0.fq12_mul(&t1, &y0);

        t0 = self.0.fq12_cyclotomic_square(&t0);
        t0 = self.0.fq12_mul(&t0, &t1);

        t0
    }

    pub fn check_pairing(
        &mut self,
        terms: &[(
            &AssignedG1Affine<C, C::Scalar>,
            &AssignedG2Affine<C, C::Scalar>,
        )],
    ) {
        let prepared_terms = terms
            .iter()
            .map(|(p, q)| (*p, self.prepare_g2(q)))
            .collect::<Vec<_>>();
        let terms = prepared_terms
            .iter()
            .map(|(p, q)| (*p, q))
            .collect::<Vec<_>>();
        let res = self.multi_miller_loop(&terms[..]);
        let res = self.final_exponentiation(&res);
        println!("{}", self.0.base_offset);
        self.0.fq12_assert_identity(&res);
    }
}
