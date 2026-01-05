//! Zero-knowledge proof support
//!
//! This module provides utilities for constructing zero-knowledge proofs,
//! including rejection sampling and Fiat-Shamir based schemes.

pub mod fiat_shamir;
pub mod rejection_sampling;

pub use fiat_shamir::{
    Challenge, Commitment, IdentificationScheme, IdentificationSchemeParameters, KeyPair,
    ProverState, PublicKey, Response, SecretKey, Signature, SignatureScheme,
};
pub use rejection_sampling::{RejectionSamplingParams, RejectionSamplingStats};
