#![cfg_attr(not(feature = "std"), no_std)]

use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Sub};

pub mod cheby;

pub trait Real: Sized + Copy + Add<Output=Self> + AddAssign + Sub<Output=Self> + Mul<Output=Self> + MulAssign + Div<Output=Self> {
    const ZERO: Self;
    const HALF: Self;
    const ONE: Self;
    const TWO: Self;
    const FOUR: Self;
    const PI: Self;

    fn from_usize(value: usize) -> Self;
}

impl Real for f32 {
    const ZERO: Self = 0.0;
    const HALF: Self = 0.5;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;
    const FOUR: Self = 4.0;
    const PI: Self = core::f32::consts::PI;

    fn from_usize(value: usize) -> Self {
        value as _
    }
}

impl Real for f64 {
    const ZERO: Self = 0.0;
    const HALF: Self = 0.5;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;
    const FOUR: Self = 4.0;
    const PI: Self = core::f64::consts::PI;

    fn from_usize(value: usize) -> Self {
        value as _
    }
}
