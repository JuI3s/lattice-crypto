//! Vector operations over Z_q

use rand::Rng;
use std::fmt;

/// A vector in Z_q^n
#[derive(Clone, Debug, PartialEq)]
pub struct Vector {
    pub coeffs: Vec<i64>,
    pub modulus: i64,
}

impl Vector {
    pub fn new(coeffs: Vec<i64>, modulus: i64) -> Self {
        let reduced: Vec<i64> = coeffs
            .into_iter()
            .map(|c| ((c % modulus) + modulus) % modulus)
            .collect();
        Vector {
            coeffs: reduced,
            modulus,
        }
    }

    pub fn zero(n: usize, modulus: i64) -> Self {
        Vector {
            coeffs: vec![0; n],
            modulus,
        }
    }

    pub fn random<R: Rng>(rng: &mut R, n: usize, bound: i64, modulus: i64) -> Self {
        let coeffs: Vec<i64> = (0..n).map(|_| rng.gen_range(-bound..=bound)).collect();
        Vector::new(coeffs, modulus)
    }

    pub fn random_ternary<R: Rng>(rng: &mut R, n: usize, modulus: i64) -> Self {
        // Ternary vector: coefficients in {-1, 0, 1}
        let coeffs: Vec<i64> = (0..n).map(|_| rng.gen_range(-1..=1)).collect();
        Vector::new(coeffs, modulus)
    }

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// ell-infinity norm (maximum absolute value of centered representatives)
    pub fn ell_inf_norm(&self) -> i64 {
        self.coeffs
            .iter()
            .map(|&c| {
                // Center the coefficient to [-q/2, q/2]
                let centered = if c > self.modulus / 2 {
                    c - self.modulus
                } else {
                    c
                };
                centered.abs()
            })
            .max()
            .unwrap_or(0)
    }

    /// Centered representative (coefficients in [-q/2, q/2])
    pub fn centered(&self) -> Vec<i64> {
        self.coeffs
            .iter()
            .map(|&c| {
                if c > self.modulus / 2 {
                    c - self.modulus
                } else {
                    c
                }
            })
            .collect()
    }

    /// Add two vectors
    pub fn add(&self, other: &Vector) -> Vector {
        assert_eq!(self.coeffs.len(), other.coeffs.len());
        assert_eq!(self.modulus, other.modulus);
        let coeffs: Vec<i64> = self
            .coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| (a + b) % self.modulus)
            .collect();
        Vector::new(coeffs, self.modulus)
    }

    /// Scalar multiplication
    pub fn scalar_mul(&self, scalar: i64) -> Vector {
        let coeffs: Vec<i64> = self.coeffs.iter().map(|&c| c * scalar).collect();
        Vector::new(coeffs, self.modulus)
    }

    /// Inner product
    pub fn inner_product(&self, other: &Vector) -> i64 {
        assert_eq!(self.coeffs.len(), other.coeffs.len());
        let sum: i64 = self
            .coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a * b)
            .sum();
        ((sum % self.modulus) + self.modulus) % self.modulus
    }
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} (mod {})", self.centered(), self.modulus)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector_operations() {
        let q = 7681;
        let v1 = Vector::new(vec![1, 2, 3, 4], q);
        let v2 = Vector::new(vec![5, 6, 7, 8], q);

        let sum = v1.add(&v2);
        assert_eq!(sum.coeffs, vec![6, 8, 10, 12]);

        let scaled = v1.scalar_mul(2);
        assert_eq!(scaled.coeffs, vec![2, 4, 6, 8]);

        let inner = v1.inner_product(&v2);
        assert_eq!(inner, 5 + 12 + 21 + 32); // 70
    }

    #[test]
    fn test_vector_linf_norm() {
        let q = 101;
        // Coefficients that wrap around
        let v = Vector::new(vec![100, 50, 1], q); // 100 mod 101 = -1 centered
        assert_eq!(v.ell_inf_norm(), 50);

        let v2 = Vector::new(vec![-10, 20, -30], q);
        assert_eq!(v2.ell_inf_norm(), 30);
    }
}
