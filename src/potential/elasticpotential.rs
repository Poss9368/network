use super::Potential;

#[derive(Debug, Clone, Copy)]
pub struct ElasticPotential {
    pub stiffness: f64,
}

impl ElasticPotential {
    pub fn new(stiffness: f64) -> Self {
        Self { stiffness }
    }
}
impl Potential for ElasticPotential {
    fn force(&self, length: f64, initial_length: f64) -> f64 {
        self.stiffness * (length - initial_length)
    }

    fn energy(&self, length: f64, initial_length: f64) -> f64 {
        0.5 * self.stiffness * (length - initial_length).powi(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elastic_energy_and_force_for_stretch_compression_and_rest() {
        let potential = ElasticPotential::new(4.0);
        for (length, force, energy) in [(1.5, 2.0, 0.5), (0.5, -2.0, 0.5), (1.0, 0.0, 0.0)] {
            assert!((potential.force(length, 1.0) - force).abs() < 1e-12);
            assert!((potential.energy(length, 1.0) - energy).abs() < 1e-12);
        }
    }

    #[test]
    fn force_matches_energy_derivative() {
        let potential = ElasticPotential::new(3.0);
        let h = 1e-6;
        for length in [0.7, 1.0, 1.4] {
            let derivative =
                (potential.energy(length + h, 1.0) - potential.energy(length - h, 1.0)) / (2.0 * h);
            assert!((potential.force(length, 1.0) - derivative).abs() < 1e-9);
        }
    }
}
