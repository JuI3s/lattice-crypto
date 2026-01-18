//! Core lattice operations and types
//!
//! This module provides fundamental data structures for working with lattices,
//! including vectors and matrices over Z_q.

pub mod matrix;
pub mod splitting;
pub mod trace;
pub mod vector;

pub use matrix::Matrix;
pub use splitting::{Challenge, ChallengeSet, SplittingParams};
pub use vector::Vector;
