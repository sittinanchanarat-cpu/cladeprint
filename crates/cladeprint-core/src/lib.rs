//! cladeprint core: near-miss filters for group-specific k-mer markers, unitigs, scoring.
//!
//! Port of the Python prototype in `kmers/`; the parity tests in `tests/python_parity.rs` check the
//! results k-mer for k-mer against its outputs. K-mer counting stays with KMC for now.

pub mod classify;
pub mod clusters;
pub mod config;
pub mod count;
pub mod edit;
pub mod groups;
pub mod hamming;
mod index;
pub mod kmc;
pub mod kmer;
pub mod npy;
pub mod pipeline;
pub mod score;
pub mod subsample;
pub mod unitig;

#[doc(hidden)]
pub mod testutil;
