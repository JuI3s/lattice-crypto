# Cyclotomic Rings

## Definition

A **cyclotomic ring** is a quotient ring of the form:

\\[ R = \mathbb{Z}[X] / \Phi_m(X) \\]

where \\( \Phi_m(X) \\) is the \\( m \\)-th cyclotomic polynomial. When working modulo a prime \\( q \\), we use:

\\[ R_q = \mathbb{Z}_q[X] / \Phi_m(X) \\]

## Power-of-Two Cyclotomics

For efficiency in lattice-based cryptography, we typically use **power-of-two cyclotomics** where \\( m = 2n \\) for some power of two \\( n \\). In this case:

\\[ \Phi_{2n}(X) = X^n + 1 \\]

This gives us the ring:

\\[ R_q = \mathbb{Z}_q[X] / (X^n + 1) \\]

### Why Power-of-Two?

1. **Efficient multiplication**: The Number Theoretic Transform (NTT) enables \\( O(n \log n) \\) polynomial multiplication
2. **Simple structure**: Elements are polynomials of degree < n with coefficients in \\( \mathbb{Z}_q \\)

## Ring Elements

An element \\( a \in R_q \\) can be written as:

\\[ a = a_0 + a_1 X + a_2 X^2 + \cdots + a_{n-1} X^{n-1} \\]

where each \\( a_i \in \mathbb{Z}_q \\).

### Norms

We define several norms on ring elements:

- **\\( \ell_\infty \\) norm**: \\( \|a\|_\infty = \max_i |a_i| \\)
- **\\( \ell_1 \\) norm**: \\( \|a\|_1 = \sum_i |a_i| \\)
- **\\( \ell_2 \\) norm**: \\( \|a\|_2 = \sqrt{\sum_i a_i^2} \\)

## Splitting Behavior

The polynomial \\( X^n + 1 \\) may factor differently modulo different primes \\( q \\):

- **Fully splits**: \\( X^n + 1 \equiv \prod_{i=1}^{n} (X - \zeta_i) \pmod{q} \\)
- **Partially splits**: Factors into \\( k \\) irreducible polynomials where \\( 1 < k < n \\)
- **Remains irreducible**: \\( X^n + 1 \\) is irreducible mod \\( q \\)

The splitting behavior depends on the multiplicative order of \\( q \\) modulo \\( 2n \\).

## Implementation

In this library, ring elements are represented by the `CyclotomicRingElement` struct:

```rust
pub struct CyclotomicRingElement {
    pub coefficients: Vec<i64>,
    pub n: usize,  // Ring dimension (X^n + 1)
}
```

For sparse polynomials (like challenges), we use a sparse representation:

```rust
pub struct Challenge {
    pub sparse_coeffs: Vec<(usize, i64)>,  // (index, coefficient) pairs
    pub n: usize,
}
```

