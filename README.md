# ternary-renormalize

Coarse-graining / renormalization group for ternary populations.

Block transformation (groups of 3→1), majority-rule coarse-graining, and tracking of quantities across scales. Inspired by statistical mechanics renormalization group methods applied to ternary-state systems.

## Features

- **5 coarse-graining rules**: Majority, SumThreshold, First, Middle, Last
- **RG flow tracking**: Automatic coarse-graining across all scales
- **Observable computation**: Fractions, entropy, magnetization, Ω at each scale
- **Fixed point detection**: Identify when RG flow converges
- **Correlation length**: Estimate from RG scale data
- **Ω stability analysis**: Check if zero-fraction is preserved under rescaling

## Usage

```rust
use ternary_renormalize::{Renormalizer, CoarseRule, Ternary, uniform_population, mixed_population};

let pop = mixed_population(243, 0.33, 0.34);
let mut renorm = Renormalizer::new(CoarseRule::Majority);
let history = renorm.rg_flow(&pop, 1);

for obs in history {
    println!("Scale {}: Ω={}, entropy={:.3}", obs.scale, obs.omega, obs.entropy);
}

// Check fixed points and stability
println!("Fixed point: {}", renorm.is_fixed_point());
println!("Ω stable: {}", renorm.omega_stable(0.01));
```

## Test Coverage

23 tests covering uniform populations, all 5 coarse-graining rules, RG flow on 81 and 243-element populations, correlation length, Ω stability, magnetization, and edge cases.

## Known Limitations

- Block size fixed at 3 (standard for ternary systems)
- Majority rule uses simple plurality (ties → Zero), not weighted voting
- No temperature-like parameter for probabilistic coarse-graining
- Triangle inequality not enforced in coarse-graining
- Fixed point detection uses exact equality within tolerance, not convergence rate
- Maximum dimension of RG flow not bounded (runs until population < min_size)

## License

MIT
