//! Rejection sampling for zero-knowledge proofs
//!
//! Implementation of the rejection sampling framework from Lyubashevsky's 2009 paper
//! "Fiat-Shamir with Aborts: Applications to Lattice and Factoring-Based Signatures"
//!
//! The key insight is that in lattice-based identification/signature schemes, the prover's
//! response z = y + sc (where y is random, s is the secret, c is the challenge) leaks
//! information about s. Rejection sampling ensures that the distribution of z is
//! independent of s by only outputting z when it falls in a "safe" range.

use crate::lattice::Vector;

/// Parameters for rejection sampling
#[derive(Clone, Debug)]
pub struct RejectionSamplingParams {
    /// Bound for the commitment randomness y: sampled from [-b, b]
    pub b: i64,
    /// Bound for the challenge coefficients
    pub challenge_bound: i64,
    /// Bound for the secret key coefficients
    pub secret_bound: i64,
    /// Dimension of vectors
    pub n: usize,
}

impl RejectionSamplingParams {
    /// Compute the rejection threshold
    /// z is accepted if ||z||_∞ <= b - challenge_bound * secret_bound * n
    pub fn rejection_bound(&self) -> i64 {
        // Conservative bound: account for worst-case ||sc||_∞
        self.b - self.challenge_bound * self.secret_bound * (self.n as i64)
    }
}

/// Statistics for rejection sampling
#[derive(Clone, Debug, Default)]
pub struct RejectionSamplingStats {
    pub total_attempts: u64,
    pub aborts: u64,
    pub successes: u64,
}

impl RejectionSamplingStats {
    pub fn abort_rate(&self) -> f64 {
        if self.total_attempts == 0 {
            0.0
        } else {
            self.aborts as f64 / self.total_attempts as f64
        }
    }
}

/// Check if a response vector should be accepted or rejected
pub fn should_accept(z: &Vector, params: &RejectionSamplingParams) -> bool {
    let rejection_bound = params.rejection_bound();
    let z_norm = z.ell_inf_norm();
    z_norm <= rejection_bound
}
