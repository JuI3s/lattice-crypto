# Trace Computation

## Definition

The **trace** of an element \\( x \in R_q \\) with respect to a Galois subgroup \\( H \\) is:

\\[ \text{Tr}_H(x) = \sum_{\sigma \in H} \sigma(x) \\]

## Naive Algorithm

The straightforward approach computes all \\( |H| \\) automorphisms and sums:

```rust
pub fn trace_naive(&self, x: &CyclotomicRingElement, h: &GaloisSubgroup) -> CyclotomicRingElement {
    let mut result = CyclotomicRingElement::zero(x.n);
    for sigma_i in h.elements() {
        let sigma_x = self.apply_automorphism(x, sigma_i);
        result = result.add(&sigma_x);
    }
    result
}
```

**Complexity**: \\( O(|H| \cdot n) \\) — one automorphism application per group element.

## Tower Optimization

When \\( H \\) admits a **tower decomposition** with tower generators \\( g_1, g_2, \ldots, g_t \\):

\\[ H = H_0 \supset H_1 \supset \cdots \supset H_t = \\{1\\} \\]

we can compute the trace iteratively:

```rust
pub fn trace_tower(&self, x: &CyclotomicRingElement, tower_gens: &[i64]) -> CyclotomicRingElement {
    let mut acc = x.clone();
    
    for &gen in tower_gens {
        // acc = acc + σ_gen(acc)
        let sigma_acc = self.apply_automorphism(&acc, gen);
        acc = acc.add(&sigma_acc);
    }
    
    acc
}
```

**Complexity**: \\( O(t \cdot n) = O(\log |H| \cdot n) \\)

## Why It Works

At each step, we're computing a partial trace:

1. Start: \\( x \\)
2. After \\( g_1 \\): \\( x + \sigma_{g_1}(x) = \text{Tr}_{H_1 \to H_0}(x) \\)
3. After \\( g_2 \\): \\( \text{Tr}_{H_2 \to H_0}(x) \\)
4. ...
5. After \\( g_t \\): \\( \text{Tr}_H(x) \\)

Each generator "doubles" the set of automorphisms being summed.

## Worked Example

Let \\( n = 8 \\), \\( H = \langle \sigma_{-1}, \sigma_5 \rangle \\) (order 4).

Tower generators: \\( [5, -1] \\)

Starting with \\( x = 1 + 2X + 3X^2 \\):

**Step 1**: Apply \\( \sigma_5 \\) and add
- \\( \sigma_5(x) = 1 + 2X^5 + 3X^{10} = 1 + 2X^5 - 3X^2 \\)
- \\( x + \sigma_5(x) = 2 + 4X^5 \\)

**Step 2**: Apply \\( \sigma_{-1} \\) and add
- Let \\( y = 2 + 4X^5 \\)
- \\( \sigma_{-1}(y) = 2 + 4X^{-5} = 2 - 4X^3 \\)
- \\( y + \sigma_{-1}(y) = 4 + 4X^5 - 4X^3 \\)

This equals \\( \sum_{\sigma \in H} \sigma(x) \\).

## Performance Comparison

| Method | Operations | For |H| = 256 |
|--------|------------|-----|
| Naive | O(\|H\| · n) | 256 automorphisms |
| Tower | O(log\|H\| · n) | 8 automorphisms |

**Speedup**: ~32× for typical parameters.

## Code

```rust
impl TraceComputer {
    pub fn compute_trace(&self, x: &CyclotomicRingElement, h: &GaloisSubgroup) -> CyclotomicRingElement {
        if let Some(tower) = h.tower_generators() {
            self.trace_tower(x, &tower)
        } else {
            self.trace_naive(x, h)
        }
    }
}
```

