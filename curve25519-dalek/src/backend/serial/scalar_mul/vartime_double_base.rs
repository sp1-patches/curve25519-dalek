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
use crate::scalar::Scalar;
use crate::traits::Identity;
use crate::window::NafLookupTable5;
use alloc::vec::Vec;
use crate::constants::ED25519_BASEPOINT_POINT;

#[cfg(not(all(target_os = "zkvm", target_vendor = "succinct")))]
/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    let a_naf = a.non_adjacent_form(5);

    #[cfg(feature = "precomputed-tables")]
    let b_naf = b.non_adjacent_form(8);
    #[cfg(not(feature = "precomputed-tables"))]
    let b_naf = b.non_adjacent_form(5);

    // Find starting index
    let mut i: usize = 255;
    for j in (0..256).rev() {
        i = j;
        if a_naf[i] != 0 || b_naf[i] != 0 {
            break;
        }
    }

    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    #[cfg(feature = "precomputed-tables")]
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    #[cfg(not(feature = "precomputed-tables"))]
    let table_B =
        &NafLookupTable5::<ProjectiveNielsPoint>::from(&constants::ED25519_BASEPOINT_POINT);

    let mut r = ProjectivePoint::identity();
    loop {
        let mut t = r.double();

        match a_naf[i].cmp(&0) {
            Ordering::Greater => t = &t.as_extended() + &table_A.select(a_naf[i] as usize),
            Ordering::Less => t = &t.as_extended() - &table_A.select(-a_naf[i] as usize),
            Ordering::Equal => {}
        }

        match b_naf[i].cmp(&0) {
            Ordering::Greater => t = &t.as_extended() + &table_B.select(b_naf[i] as usize),
            Ordering::Less => t = &t.as_extended() - &table_B.select(-b_naf[i] as usize),
            Ordering::Equal => {}
        }

        r = t.as_projective();

        if i == 0 {
            break;
        }
        i -= 1;
    }

    r.as_extended()
}

#[cfg(all(target_os = "zkvm", target_vendor = "succinct"))]
use sp1_lib::{ed25519::Ed25519AffinePoint, utils::AffinePoint};
#[cfg(all(target_os = "zkvm", target_vendor = "succinct"))]
/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
///
/// Accelerated with SP1's EdAdd syscall.
#[allow(non_snake_case)]
pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    let A: Ed25519AffinePoint = (*A).into();

    let a_bits = a.bits_le();
    let a_bits = a_bits.iter().map(|bit| *bit == 1).collect::<Vec<bool>>();
    let b_bits = b.bits_le();
    let b_bits = b_bits.iter().map(|bit| *bit == 1).collect::<Vec<bool>>();

    // Note: The base point is the identity point.
    let res = AffinePoint::multi_scalar_multiplication(
        &a_bits,
        A,
        &b_bits,
        ED25519_BASEPOINT_POINT.into(),
    )
    .unwrap();
    res.into()
}
