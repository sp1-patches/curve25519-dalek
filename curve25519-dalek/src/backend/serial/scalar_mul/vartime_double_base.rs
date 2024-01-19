// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>
#![allow(non_snake_case)]

use core::cmp::Ordering;

use crate::backend::serial::curve_models::{ProjectiveNielsPoint, ProjectivePoint};
use crate::constants;
use crate::edwards::EdwardsPoint;
use crate::field::FieldElement;
use crate::scalar::Scalar;
use crate::traits::Identity;
use crate::window::NafLookupTable5;
use alloc::vec::Vec;

struct AffinePoint {
    base_ptr: *mut u8,
}

impl Clone for AffinePoint {
    fn clone(&self) -> Self {
        let mut serialized = [0u8; 64];
        unsafe {
            core::ptr::copy_nonoverlapping(self.base_ptr, serialized.as_mut_ptr(), 64);
        }
        AffinePoint {
            base_ptr: serialized.as_mut_ptr(),
        }
    }
}

struct AffinePointNormal {
    pub x: FieldElement,
    pub y: FieldElement,
}

impl AffinePointNormal {
    pub fn add_assign(&self, other: &AffinePointNormal) -> AffinePointNormal {
        let x1 = self.x;
        let x2 = other.x;
        let y1 = self.y;
        let y2 = other.y;

        let x3_numerator = &(&x1 * &y2) + &(&x2 * &y1);
        let y3_numerator = &(&y1 * &y2) + &(&x1 * &x2);
        let f = &(&x1 * &x2) * &(&y1 * &y2);
        let d_mul_f = &f * &constants::EDWARDS_D;

        let x3 = &x3_numerator * &(&FieldElement::ONE + &d_mul_f).invert();
        let y3 = &y3_numerator * &(&FieldElement::ONE - &d_mul_f).invert();

        AffinePointNormal { x: x3, y: y3 }
    }

    pub fn from_edwards(p: EdwardsPoint) -> Self {
        let z_inv = p.Z.invert(); // Assuming FieldElement has an invert method
        let x_affine = &p.X * &z_inv;
        let y_affine = &p.Y * &z_inv;

        AffinePointNormal {
            x: x_affine,
            y: y_affine,
        }
    }

    pub fn to_edwards(&self) -> EdwardsPoint {
        EdwardsPoint {
            X: self.x,
            Y: self.y,
            Z: FieldElement::ONE,
            T: &self.x * &self.y,
        }
    }

    pub fn double(&mut self) {
        // ????
        let x1 = self.x;
        let y1 = self.y;

        let x3_numerator = &(&x1 * &y1) + &(&x1 * &y1);
        let y3_numerator = &(&y1 * &y1) + &(&x1 * &x1);
        let f = &(&x1 * &x1) * &(&y1 * &y1);
        let d_mul_f = &f * &constants::EDWARDS_D;

        let x3 = &x3_numerator * &(&FieldElement::ONE + &d_mul_f).invert();
        let y3 = &y3_numerator * &(&FieldElement::ONE - &d_mul_f).invert();

        self.x = x3;
        self.y = y3;
    }
}

impl AffinePoint {
    pub fn from_edwards(p: EdwardsPoint) -> Self {
        let z_inv = p.Z.invert(); // Assuming FieldElement has an invert method
        let x_affine = &p.X * &z_inv;
        let y_affine = &p.Y * &z_inv;

        let mut serialized = [0u8; 64];
        serialized[..32].copy_from_slice(&x_affine.as_bytes());
        serialized[32..].copy_from_slice(&y_affine.as_bytes());

        AffinePoint {
            base_ptr: serialized.as_mut_ptr(),
        }
    }

    pub fn set_normal(&mut self, p: AffinePointNormal) {
        let bytes = p.x.as_bytes();
        unsafe {
            core::ptr::copy_nonoverlapping(bytes.as_ptr(), self.base_ptr, 32);
        }
        let bytes = p.y.as_bytes();
        unsafe {
            core::ptr::copy_nonoverlapping(bytes.as_ptr(), self.base_ptr.add(32), 32);
        }
    }

    pub fn from_normal(p: AffinePointNormal) -> Self {
        let mut serialized = [0u8; 64];
        serialized[..32].copy_from_slice(&p.x.as_bytes());
        serialized[32..].copy_from_slice(&p.y.as_bytes());

        AffinePoint {
            base_ptr: serialized.as_mut_ptr(),
        }
    }

    pub fn to_normal(&self) -> AffinePointNormal {
        let mut bytes: [u8; 32] = [0u8; 32];
        unsafe {
            core::ptr::copy_nonoverlapping(self.base_ptr, bytes.as_mut_ptr(), 32);
        }
        let x = FieldElement::from_bytes(&bytes);
        unsafe {
            core::ptr::copy_nonoverlapping(self.base_ptr.add(32), bytes.as_mut_ptr(), 32);
        }
        let y = FieldElement::from_bytes(&bytes);

        AffinePointNormal { x, y }
    }

    // TODO implement Add with another AffinePoint
    pub fn add_assign(&mut self, other: &AffinePoint) {
        let a = self.to_normal();
        let b = other.to_normal();
        let c = a.add_assign(&b);
        self.set_normal(c);
        // ecall to ed_add
    }

    pub fn double(&mut self) {
        let a = self.to_normal();
        let b = a.add_assign(&a);
        self.set_normal(b);
        // ecall to ed_double
    }

    pub fn to_edwards(&self) -> EdwardsPoint {
        let mut bytes: [u8; 32] = [0u8; 32];
        unsafe {
            core::ptr::copy_nonoverlapping(self.base_ptr, bytes.as_mut_ptr(), 32);
        }
        let x = FieldElement::from_bytes(&bytes);
        unsafe {
            core::ptr::copy_nonoverlapping(self.base_ptr.add(32), bytes.as_mut_ptr(), 32);
        }
        let y = FieldElement::from_bytes(&bytes);

        EdwardsPoint {
            X: x,
            Y: y,
            Z: FieldElement::ONE,
            T: &x * &y,
        }
    }
}

pub fn ecall_mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    let mut a_decomp: Vec<bool> = a.bits_le().collect();
    let mut b_decomp: Vec<bool> = b.bits_le().collect();
    if a_decomp.len() < b_decomp.len() {
        a_decomp.resize(b_decomp.len(), false);
    } else {
        b_decomp.resize(a_decomp.len(), false);
    }
    let max_len = a_decomp.len();
    let mut temp_A = AffinePointNormal::from_edwards(*A);
    let mut temp_B = AffinePointNormal::from_edwards(constants::ED25519_BASEPOINT_POINT);

    let mut res = AffinePointNormal::from_edwards(EdwardsPoint::identity());

    for bit in 0..max_len {
        if a_decomp[bit] == true {
            res.add_assign(&temp_A);
        }
        if b_decomp[bit] == true {
            res.add_assign(&temp_B);
        }
        temp_A.double();
        temp_B.double();
    }
    res.to_edwards()
}

// pub fn ecall_mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
//     let mut a_decomp: Vec<bool> = a.bits_le().collect();
//     let mut b_decomp: Vec<bool> = b.bits_le().collect();
//     if a_decomp.len() < b_decomp.len() {
//         a_decomp.resize(b_decomp.len(), false);
//     } else {
//         b_decomp.resize(a_decomp.len(), false);
//     }
//     let max_len = a_decomp.len();
//     let mut temp_A = AffinePoint::from_edwards(*A);
//     let mut temp_B = AffinePoint::from_edwards(constants::ED25519_BASEPOINT_POINT);

//     let mut res = AffinePoint::from_edwards(EdwardsPoint::identity());

//     for bit in 0..max_len {
//         if a_decomp[bit] == true {
//             res.add_assign(&temp_A);
//         }
//         if b_decomp[bit] == true {
//             res.add_assign(&temp_B);
//         }
//         temp_A.double();
//         temp_B.double();
//     }
//     return res.to_edwards();
// }

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    ecall_mul(a, A, b)
    // let a_naf = a.non_adjacent_form(5);

    // #[cfg(feature = "precomputed-tables")]
    // let b_naf = b.non_adjacent_form(8);
    // #[cfg(not(feature = "precomputed-tables"))]
    // let b_naf = b.non_adjacent_form(5);

    // // Find starting index
    // let mut i: usize = 255;
    // for j in (0..256).rev() {
    //     i = j;
    //     if a_naf[i] != 0 || b_naf[i] != 0 {
    //         break;
    //     }
    // }

    // let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    // #[cfg(feature = "precomputed-tables")]
    // let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    // #[cfg(not(feature = "precomputed-tables"))]
    // let table_B =
    //     &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_POINT);

    // let mut r = ProjectivePoint::identity();
    // loop {
    //     let mut t = r.double();

    //     match a_naf[i].cmp(&0) {
    //         Ordering::Greater => t = &t.as_extended() + &table_A.select(a_naf[i] as usize),
    //         Ordering::Less => t = &t.as_extended() - &table_A.select(-a_naf[i] as usize),
    //         Ordering::Equal => {}
    //     }

    //     match b_naf[i].cmp(&0) {
    //         Ordering::Greater => t = &t.as_extended() + &table_B.select(b_naf[i] as usize),
    //         Ordering::Less => t = &t.as_extended() - &table_B.select(-b_naf[i] as usize),
    //         Ordering::Equal => {}
    //     }

    //     r = t.as_projective();

    //     if i == 0 {
    //         break;
    //     }
    //     i -= 1;
    // }

    // r.as_extended()
}
