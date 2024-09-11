use approximate::cheby::Chebyshev;
use std::f32::consts::PI;
use std::mem;

fn main() {
    let min_bound = 0.0;
    let max_bound = 1.0;
    let f = |x: f32| x.cos();
    let cheby = Chebyshev::<6, _>::fit(min_bound, max_bound, f, |x| x.cos());
    let g = |x: f32| { cheby.eval(x) };
    let mut total_err: u64 = 0;
    let mut max_err: u32 = 0;
    for err in iter_f32(min_bound, max_bound).map(|x| float_distance(f(x), g(x))) {
        total_err += err as u64;
        max_err = max_err.max(err);
    }
    dbg!(max_err);
    dbg!(total_err);
    dbg!(num_floats(min_bound, max_bound) as f64);
    dbg!(total_err as f64 / num_floats(min_bound, max_bound) as f64);
    dbg!(cheby);
}

// cargo asm -p approximate --example=floats --simplify --intel asm_eval
#[no_mangle]
#[inline(never)]
pub fn asm_eval(x: f32) -> f32 {
    const CHEBY: Chebyshev<8, f32> = Chebyshev::from_coeffs(0.0, 1.2732395, [
        -0.2498861,
        -1.0473543,
        0.27981228,
        0.052837715,
        -0.036263496,
        0.0017184913,
        -0.00045645982,
        -6.323308e-5,
    ]);
    CHEBY.eval(x)
}

/// Distance between two floats.
/// e.g. returns 0 if a == b, 1 if a and b are consecutive ...
fn float_distance(a: f32, b: f32) -> u32 {
    debug_assert_eq!(a.is_sign_positive(), b.is_sign_positive(), "min and max must be of the same sign");
    debug_assert!(!a.is_nan() && !a.is_infinite());
    debug_assert!(!b.is_nan() && !b.is_infinite());
    let max = unsafe { mem::transmute::<_, u32>(b) };
    let min = unsafe { mem::transmute::<_, u32>(a) };
    max.abs_diff(min)
}

fn num_floats(a: f32, b: f32) -> u32 {
    float_distance(a, b) + 1
}

struct Iterf32 {
    float: u32,
    max: f32,
}

impl Iterf32 {
    fn read_float(&self) -> f32 {
        unsafe { mem::transmute::<_, f32>(self.float) }
    }
}

impl Iterator for Iterf32 {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        let value = self.read_float();
        if value <= self.max {
            self.float += 1;
            Some(value)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let num = num_floats(self.read_float(), self.max) as usize;
        (num, Some(num))
    }
}

impl ExactSizeIterator for Iterf32 {}

fn iter_f32(min: f32, max: f32) -> impl ExactSizeIterator<Item=f32> {
    Iterf32 { float: unsafe { mem::transmute::<_, _>(min) }, max }
}
