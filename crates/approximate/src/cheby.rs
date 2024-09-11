use crate::Real;

pub fn dct_coeff<T: Real>(cos: impl Fn(T) -> T, N: usize, outer: usize, inner: usize) -> T {
    cos(T::PI * T::from_usize(outer) * (T::from_usize(inner) + T::HALF) / T::from_usize(N))
}

fn gauss_chebyshev_zeros<T: Real>(cos: impl Fn(T) -> T, N: usize, k: usize) -> T {
    cos(T::PI * (T::from_usize(k) + T::HALF) / T::from_usize(N))
}

// jeez, mathematicians
const fn kronecker_delta<T: Real>(i: usize, j: usize) -> T {
    if i == j {
        T::ONE
    } else {
        T::ZERO
    }
}

#[cfg_attr(feature = "std", derive(Debug))]
pub struct Chebyshev<const N: usize, T: Real> {
    x_min: T,
    x_scale: T,
    coeffs: [T; N],
}

impl<const N: usize, T: Real> Chebyshev<N, T> {
    pub const fn from_coeffs(x_min: T, x_scale: T, coeffs: [T; N]) -> Self {
        Self {
            x_min,
            x_scale,
            coeffs,
        }
    }

    /// N evaluations of f and N^2 evaluation of cos
    pub fn fit(x_min: T, x_max: T, f: impl Fn(T) -> T, cos: impl Fn(T) -> T + Copy) -> Self {
        // Instead of using the inner product method which could require us to compute integrals,
        // We can rely on the following property : Tn(cos(x)) = cos(nx) and by application of the discrete orthogonality condition
        // With xk the Gauss-Chebyshev zeros of Tn, Tn(xk) = Tn(cos(pi(k + 1/2)/N) = cos(n pi(k + 1/2)/N) = dct coefficient
        // https://en.wikipedia.org/wiki/Chebyshev_polynomials#As_a_basis_set
        // I guess that if N is large enough, we could use a fast cosine transform to compute the coefficients,
        // but it's not like we are calling 'fit' all over the place.
        let x_range = x_max - x_min;
        let mut values = [T::ZERO; N];
        for k in 0..N {
            // remap the zeros from [-1,1] to [x_min,x_max]
            values[k] = f((gauss_chebyshev_zeros(cos, N, k) + T::ONE) * x_range * T::HALF + x_min);
        }
        // let cos = cos;
        let mut coeffs = [T::ZERO; N];
        for n in 0..N {
            for k in 0..N {
                coeffs[n] += values[k] * dct_coeff(cos, N, n, k);
            }
            coeffs[n] *= (T::TWO - kronecker_delta(0, n)) / T::from_usize(N);
        }

        Self { x_min, x_scale: T::FOUR / x_range, coeffs }
    }

    pub fn fit_match(x_min: T, x_max: T, f: impl Fn(T) -> T, cos: impl Fn(T) -> T + Copy) -> Self {
        // Instead of using the inner product method which could require us to compute integrals,
        // We can rely on the following property : Tn(cos(x)) = cos(nx) and by application of the discrete orthogonality condition
        // With xk the Gauss-Chebyshev zeros of Tn, Tn(xk) = Tn(cos(pi(k + 1/2)/N) = cos(n pi(k + 1/2)/N) = dct coefficient
        // https://en.wikipedia.org/wiki/Chebyshev_polynomials#As_a_basis_set
        // I guess that if N is large enough, we could use a fast cosine transform to compute the coefficients,
        // but it's not like we are calling 'fit' all over the place.
        let x_range = x_max - x_min;
        let mut approx = Self::from_coeffs(x_min, T::FOUR / x_range, [T::ZERO; N]);
        let mut values = [T::ZERO; N];
        for k in 0..N {
            // remap the zeros from [-1,1] to [x_min,x_max]
            values[k] = f((gauss_chebyshev_zeros(cos, N, k) + T::ONE) * x_range * T::HALF + x_min);
        }
        // let cos = cos;
        for n in 0..N {
            for k in 0..N {
                approx.coeffs[n] += values[k] * dct_coeff(cos, N, n, k);
            }
            approx.coeffs[n] *= (T::TWO - kronecker_delta(0, n)) / T::from_usize(N);
        }

        // Add a linear term a + bx that offsets the left and right
        // ends to the desired values
        let (x_min_offs, x_max_offs) = (f(x_min) - approx.eval(x_min), f(x_max) - approx.eval(x_max));
        let a = T::HALF * (x_max_offs + x_min_offs);
        let b = T::HALF * (x_max_offs - x_min_offs);
        approx.coeffs[0] += a; // multiplied by 2.0 * 0.5 = 1 due to c0 pre-bake multiply above
        approx.coeffs[1] += b;

        approx
    }

    pub fn eval(&self, x: T) -> T {
        // Chebyshev polynomials are defined recursively which is bad for computations
        // We use the Clenshaw algorithm to easily evaluate our approximation.
        // https://en.wikipedia.org/wiki/Clenshaw_algorithm

        // - Pre-multiply once by 2 so we don't have to do it in the loop
        // - We also need to remap x from [x_min,x_max] to [-1,1]
        // We remap x to [-2,2] to account for both statements
        let x = (x - self.x_min) * self.x_scale - T::TWO;
        let mut bk = T::ZERO;
        let mut bk1 = T::ZERO;
        let mut bk2 = T::ZERO;

        let mut k = N - 1;
        while k > 0 {
            bk = self.coeffs[k] + x * bk1 - bk2;
            bk2 = bk1;
            bk1 = bk;
            k -= 1;
        }

        // Multiply by 0.5 to invert the pre-multiplication
        self.coeffs[0] + T::HALF * x * bk - bk2
    }

    pub fn eval_trunc(&self, x: T, trunc: usize) -> T {
        // Chebyshev polynomials are defined recursively which is bad for computations
        // We use the Clenshaw algorithm to easily evaluate our approximation.
        // https://en.wikipedia.org/wiki/Clenshaw_algorithm

        // - Pre-multiply once by 2 so we don't have to do it in the loop
        // - We also need to remap x from [x_min,x_max] to [-1,1]
        // We remap x to [-2,2] to account for both statements
        let x = (x - self.x_min) * self.x_scale - T::TWO;
        let mut bk = T::ZERO;
        let mut bk1 = T::ZERO;
        let mut bk2 = T::ZERO;

        let mut k = (N - 1).min(trunc);
        while k > 0 {
            bk = self.coeffs[k] + x * bk1 - bk2;
            bk2 = bk1;
            bk1 = bk;
            k -= 1;
        }

        // Multiply by 0.5 to invert the pre-multiplication
        self.coeffs[0] + T::HALF * x * bk - bk2
    }
}

// const float arithmetic
#[rustversion::since(1.82)]
impl<const N: usize> Chebyshev<N, f32> {
    pub const fn const_eval(&self, x: f32) -> f32 {
        // Chebyshev polynomials are defined recursively which is bad for computations
        // We use the Clenshaw algorithm to easily evaluate our approximation.
        // https://en.wikipedia.org/wiki/Clenshaw_algorithm

        // - Pre-multiply once by 2 so we don't have to do it in the loop
        // - We also need to remap x from [x_min,x_max] to [-1,1]
        // We remap x to [-2,2] to account for both statements
        let x = (x - self.x_min) * self.x_scale - 2.0;
        let mut bk = 0.0;
        let mut bk1 = 0.0;
        let mut bk2 = 0.0;

        let mut k = N - 1;
        while k > 0 {
            bk = self.coeffs[k] + x * bk1 - bk2;
            bk2 = bk1;
            bk1 = bk;
            k -= 1;
        }

        // Multiply by 0.5 to invert the pre-multiplication
        self.coeffs[0] + 0.5 * x * bk - bk2
    }
}

#[rustversion::since(1.82)]
impl<const N: usize> Chebyshev<N, f64> {
    pub const fn const_eval(&self, x: f64) -> f64 {
        // Chebyshev polynomials are defined recursively which is bad for computations
        // We use the Clenshaw algorithm to easily evaluate our approximation.
        // https://en.wikipedia.org/wiki/Clenshaw_algorithm

        // - Pre-multiply once by 2 so we don't have to do it in the loop
        // - We also need to remap x from [x_min,x_max] to [-1,1]
        // We remap x to [-2,2] to account for both statements
        let x = (x - self.x_min) * self.x_scale - 2.0;
        let mut bk = 0.0;
        let mut bk1 = 0.0;
        let mut bk2 = 0.0;

        let mut k = N - 1;
        while k > 0 {
            bk = self.coeffs[k] + x * bk1 - bk2;
            bk2 = bk1;
            bk1 = bk;
            k -= 1;
        }

        // Multiply by 0.5 to invert the pre-multiplication
        self.coeffs[0] + 0.5 * x * bk - bk2
    }
}
