#![forbid(unsafe_code)]

//! # ternary-renormalize
//!
//! Coarse-graining / renormalization group for ternary populations.
//!
//! Block transformation (groups of 3→1), majority-rule coarse-graining,
//! tracking quantities across scales.

/// A ternary state: -1, 0, +1
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Ternary {
    Neg,
    Zero,
    Pos,
}

impl Ternary {
    pub fn value(&self) -> i8 {
        match self {
            Ternary::Neg => -1,
            Ternary::Zero => 0,
            Ternary::Pos => 1,
        }
    }

    pub fn from_value(v: i8) -> Self {
        if v < 0 { Ternary::Neg } else if v > 0 { Ternary::Pos } else { Ternary::Zero }
    }

    pub fn to_index(&self) -> usize {
        match self {
            Ternary::Neg => 0,
            Ternary::Zero => 1,
            Ternary::Pos => 2,
        }
    }
}

/// Coarse-graining rule
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoarseRule {
    /// Majority vote among 3 sites (ties → Zero)
    Majority,
    /// Sum and threshold: sum > 0 → Pos, sum < 0 → Neg, else Zero
    SumThreshold,
    /// First element of block
    First,
    /// Middle element of block
    Middle,
    /// Last element of block
    Last,
}

/// Observable quantities tracked across RG scales
#[derive(Debug, Clone)]
pub struct RGObservable {
    /// Fraction of Neg states
    pub frac_neg: f64,
    /// Fraction of Zero states
    pub frac_zero: f64,
    /// Fraction of Pos states
    pub frac_pos: f64,
    /// Total population size at this scale
    pub population_size: usize,
    /// Scale level (0 = original, 1 = first coarse-graining, etc.)
    pub scale: usize,
    /// Entropy at this scale
    pub entropy: f64,
    /// Magnetization (mean value)
    pub magnetization: f64,
    /// Omega parameter: fraction of zeros (Ω)
    pub omega: f64,
}

impl RGObservable {
    pub fn from_population(pop: &[Ternary], scale: usize) -> Self {
        let n = pop.len() as f64;
        let neg = pop.iter().filter(|&&s| s == Ternary::Neg).count() as f64;
        let zero = pop.iter().filter(|&&s| s == Ternary::Zero).count() as f64;
        let pos = pop.iter().filter(|&&s| s == Ternary::Pos).count() as f64;

        let frac_neg = neg / n;
        let frac_zero = zero / n;
        let frac_pos = pos / n;

        let probs = [frac_neg, frac_zero, frac_pos];
        let entropy = -probs.iter()
            .filter(|&&p| p > 0.0)
            .map(|&p| p * p.log2())
            .sum::<f64>();

        let mag: f64 = pop.iter().map(|s| s.value() as f64).sum::<f64>() / n;

        RGObservable {
            frac_neg,
            frac_zero,
            frac_pos,
            population_size: pop.len(),
            scale,
            entropy,
            magnetization: mag,
            omega: frac_zero,
        }
    }
}

/// Renormalization group transformation engine
pub struct Renormalizer {
    rule: CoarseRule,
    scale_history: Vec<RGObservable>,
    population_history: Vec<Vec<Ternary>>,
}

impl Renormalizer {
    pub fn new(rule: CoarseRule) -> Self {
        Self {
            rule,
            scale_history: Vec::new(),
            population_history: Vec::new(),
        }
    }

    /// Perform a single RG step (block transformation groups of 3 → 1)
    pub fn coarse_grain(&self, pop: &[Ternary]) -> Vec<Ternary> {
        let mut result = Vec::new();
        for chunk in pop.chunks(3) {
            let mapped = match self.rule {
                CoarseRule::Majority => self.majority(chunk),
                CoarseRule::SumThreshold => {
                    let sum: i32 = chunk.iter().map(|s| s.value() as i32).sum();
                    Ternary::from_value(if sum > 0 { 1 } else if sum < 0 { -1 } else { 0 })
                }
                CoarseRule::First => chunk[0],
                CoarseRule::Middle => chunk.get(chunk.len() / 2).copied().unwrap_or(Ternary::Zero),
                CoarseRule::Last => *chunk.last().unwrap_or(&Ternary::Zero),
            };
            result.push(mapped);
        }
        result
    }

    fn majority(&self, chunk: &[Ternary]) -> Ternary {
        let neg = chunk.iter().filter(|&&s| s == Ternary::Neg).count();
        let zero = chunk.iter().filter(|&&s| s == Ternary::Zero).count();
        let pos = chunk.iter().filter(|&&s| s == Ternary::Pos).count();
        if neg > zero && neg > pos { Ternary::Neg }
        else if pos > neg && pos > zero { Ternary::Pos }
        else { Ternary::Zero }
    }

    /// Run full RG flow, coarse-graining until population < min_size
    pub fn rg_flow(&mut self, initial: &[Ternary], min_size: usize) -> &[RGObservable] {
        self.scale_history.clear();
        self.population_history.clear();

        let mut current = initial.to_vec();
        let mut scale = 0;

        loop {
            let obs = RGObservable::from_population(&current, scale);
            self.scale_history.push(obs);
            self.population_history.push(current.clone());

            if current.len() < min_size || current.len() < 3 {
                break;
            }

            current = self.coarse_grain(&current);
            scale += 1;
        }

        &self.scale_history
    }

    /// Get the scale history
    pub fn scale_history(&self) -> &[RGObservable] {
        &self.scale_history
    }

    /// Get population at a given scale
    pub fn population_at_scale(&self, scale: usize) -> Option<&[Ternary]> {
        self.population_history.get(scale).map(|v| v.as_slice())
    }

    /// Check if the RG flow has reached a fixed point
    pub fn is_fixed_point(&self) -> bool {
        if self.scale_history.len() < 2 {
            return false;
        }
        let last = &self.scale_history[self.scale_history.len() - 1];
        let prev = &self.scale_history[self.scale_history.len() - 2];
        (last.frac_neg - prev.frac_neg).abs() < 1e-10
            && (last.frac_zero - prev.frac_zero).abs() < 1e-10
            && (last.frac_pos - prev.frac_pos).abs() < 1e-10
    }

    /// Compute correlation length from RG data
    pub fn correlation_length(&self) -> f64 {
        if self.scale_history.len() < 2 {
            return 0.0;
        }
        // At each step, length scale multiplies by 3 (block size)
        // Correlation length ξ ≈ 3^n where n is the scale at which
        // the system reaches its fixed point
        let n = self.scale_history.len() - 1;
        3.0_f64.powi(n as i32)
    }

    /// Check if Omega (fraction of zeros) is stable under rescaling
    pub fn omega_stable(&self, tolerance: f64) -> bool {
        if self.scale_history.len() < 2 {
            return true;
        }
        let first_omega = self.scale_history[0].omega;
        self.scale_history.iter().all(|obs| (obs.omega - first_omega).abs() < tolerance)
    }
}

/// Generate a uniform ternary population
pub fn uniform_population(n: usize, state: Ternary) -> Vec<Ternary> {
    vec![state; n]
}

/// Generate a random-ish ternary population with given fractions
pub fn mixed_population(n: usize, frac_neg: f64, frac_zero: f64) -> Vec<Ternary> {
    let neg_count = (n as f64 * frac_neg) as usize;
    let zero_count = (n as f64 * frac_zero) as usize;
    let mut pop = Vec::with_capacity(n);
    for _ in 0..neg_count { pop.push(Ternary::Neg); }
    for _ in 0..zero_count { pop.push(Ternary::Zero); }
    for _ in neg_count + zero_count..n { pop.push(Ternary::Pos); }
    pop
}

/// Compute Omega (fraction of zeros) for a population
pub fn omega(pop: &[Ternary]) -> f64 {
    pop.iter().filter(|&&s| s == Ternary::Zero).count() as f64 / pop.len().max(1) as f64
}

/// Compute magnetization
pub fn magnetization(pop: &[Ternary]) -> f64 {
    pop.iter().map(|s| s.value() as f64).sum::<f64>() / pop.len().max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_all_pos() {
        let pop = uniform_population(27, Ternary::Pos);
        assert!(pop.iter().all(|&s| s == Ternary::Pos));
    }

    #[test]
    fn test_uniform_all_neg() {
        let pop = uniform_population(27, Ternary::Neg);
        assert!(pop.iter().all(|&s| s == Ternary::Neg));
    }

    #[test]
    fn test_uniform_all_zero() {
        let pop = uniform_population(27, Ternary::Zero);
        assert!(pop.iter().all(|&s| s == Ternary::Zero));
    }

    #[test]
    fn test_majority_rule_all_pos() {
        let pop = uniform_population(9, Ternary::Pos);
        let renorm = Renormalizer::new(CoarseRule::Majority);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result.len(), 3);
        assert!(result.iter().all(|&s| s == Ternary::Pos));
    }

    #[test]
    fn test_majority_rule_all_neg() {
        let pop = uniform_population(9, Ternary::Neg);
        let renorm = Renormalizer::new(CoarseRule::Majority);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result.len(), 3);
        assert!(result.iter().all(|&s| s == Ternary::Neg));
    }

    #[test]
    fn test_majority_rule_all_zero() {
        let pop = uniform_population(9, Ternary::Zero);
        let renorm = Renormalizer::new(CoarseRule::Majority);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result.len(), 3);
        assert!(result.iter().all(|&s| s == Ternary::Zero));
    }

    #[test]
    fn test_majority_mixed() {
        // 2 pos, 1 neg → pos
        let pop = vec![Ternary::Pos, Ternary::Pos, Ternary::Neg];
        let renorm = Renormalizer::new(CoarseRule::Majority);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result[0], Ternary::Pos);
    }

    #[test]
    fn test_majority_tie_goes_zero() {
        // 1 pos, 1 neg, 1 zero → zero (tie)
        let pop = vec![Ternary::Pos, Ternary::Neg, Ternary::Zero];
        let renorm = Renormalizer::new(CoarseRule::Majority);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result[0], Ternary::Zero);
    }

    #[test]
    fn test_sum_threshold() {
        let pop = vec![Ternary::Pos, Ternary::Pos, Ternary::Pos]; // sum=3
        let renorm = Renormalizer::new(CoarseRule::SumThreshold);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result[0], Ternary::Pos);
    }

    #[test]
    fn test_sum_threshold_neg() {
        let pop = vec![Ternary::Neg, Ternary::Neg, Ternary::Pos]; // sum=-1
        let renorm = Renormalizer::new(CoarseRule::SumThreshold);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result[0], Ternary::Neg);
    }

    #[test]
    fn test_first_rule() {
        let pop = vec![Ternary::Neg, Ternary::Pos, Ternary::Zero, Ternary::Pos, Ternary::Neg, Ternary::Zero];
        let renorm = Renormalizer::new(CoarseRule::First);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result, vec![Ternary::Neg, Ternary::Pos]);
    }

    #[test]
    fn test_rg_flow_uniform_pos() {
        let pop = uniform_population(81, Ternary::Pos);
        let mut renorm = Renormalizer::new(CoarseRule::Majority);
        renorm.rg_flow(&pop, 1);
        let history = renorm.scale_history();
        // 81 → 27 → 9 → 3 → 1 = 5 scales
        assert_eq!(history.len(), 5);
        assert!(history.iter().all(|obs| obs.frac_pos == 1.0));
    }

    #[test]
    fn test_rg_flow_uniform_is_fixed_point() {
        let pop = uniform_population(81, Ternary::Zero);
        let mut renorm = Renormalizer::new(CoarseRule::Majority);
        renorm.rg_flow(&pop, 1);
        assert!(renorm.is_fixed_point());
    }

    #[test]
    fn test_rg_flow_random_converges() {
        let pop = mixed_population(243, 0.33, 0.34);
        let mut renorm = Renormalizer::new(CoarseRule::Majority);
        renorm.rg_flow(&pop, 1);
        // Should have multiple scales
        assert!(renorm.scale_history().len() > 2);
    }

    #[test]
    fn test_correlation_length() {
        let pop = uniform_population(81, Ternary::Pos);
        let mut renorm = Renormalizer::new(CoarseRule::Majority);
        renorm.rg_flow(&pop, 1);
        let xi = renorm.correlation_length();
        assert!(xi > 0.0);
    }

    #[test]
    fn test_omega_uniform() {
        let pop = uniform_population(27, Ternary::Zero);
        assert!((omega(&pop) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_omega_no_zeros() {
        let pop = uniform_population(27, Ternary::Pos);
        assert!((omega(&pop) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_omega_stable_under_rescaling() {
        // All-zeros should have stable omega=1 under any rescaling
        let pop = uniform_population(243, Ternary::Zero);
        let mut renorm = Renormalizer::new(CoarseRule::Majority);
        renorm.rg_flow(&pop, 1);
        assert!(renorm.omega_stable(0.01));
    }

    #[test]
    fn test_magnetization() {
        let pop = vec![Ternary::Pos, Ternary::Pos, Ternary::Neg];
        let m = magnetization(&pop);
        assert!((m - 1.0/3.0).abs() < 1e-10, "Expected 1/3, got {}", m);
    }

    #[test]
    fn test_population_at_scale() {
        let pop = uniform_population(27, Ternary::Pos);
        let mut renorm = Renormalizer::new(CoarseRule::Majority);
        renorm.rg_flow(&pop, 1);
        let s0 = renorm.population_at_scale(0).unwrap();
        assert_eq!(s0.len(), 27);
        let s1 = renorm.population_at_scale(1).unwrap();
        assert_eq!(s1.len(), 9);
    }

    #[test]
    fn test_middle_rule() {
        let pop = vec![Ternary::Neg, Ternary::Zero, Ternary::Pos];
        let renorm = Renormalizer::new(CoarseRule::Middle);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result[0], Ternary::Zero);
    }

    #[test]
    fn test_last_rule() {
        let pop = vec![Ternary::Neg, Ternary::Zero, Ternary::Pos];
        let renorm = Renormalizer::new(CoarseRule::Last);
        let result = renorm.coarse_grain(&pop);
        assert_eq!(result[0], Ternary::Pos);
    }

    #[test]
    fn test_non_multiple_of_three() {
        let pop = vec![Ternary::Pos, Ternary::Pos]; // only 2 elements
        let renorm = Renormalizer::new(CoarseRule::Majority);
        let result = renorm.coarse_grain(&pop);
        // chunk of 2 → majority of 2 pos = pos
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], Ternary::Pos);
    }
}
