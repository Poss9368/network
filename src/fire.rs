//! FIRE con velocity Verlet, masas unitarias y caja periódica fija.
//! Referencia: https://doi.org/10.1103/PhysRevLett.97.170201

use crate::{network::Network, potential::Potential};
use std::io;

#[derive(Debug, Clone, Copy)]
pub struct Fire {
    pub force_tolerance: f64,
    pub max_steps: usize,
    pub initial_dt: f64,
    pub max_dt: f64,
}

impl Default for Fire {
    fn default() -> Self {
        Self {
            force_tolerance: 1e-12,
            max_steps: 1_000_000,
            initial_dt: 0.01,
            max_dt: 0.1,
        }
    }
}

#[derive(Debug)]
pub struct FireResult {
    pub converged: bool,
    pub steps: usize,
    pub energy: f64,
    /// Máxima norma de la fuerza total sobre una partícula.
    pub max_force: f64,
}

impl Fire {
    /// Relaja posiciones; conserva caja, contactos, potencial y longitudes de reposo.
    /// Al agotar max_steps devuelve converged=false, no un éxito de minimización.
    pub fn minimize<P: Potential>(&self, network: &mut Network<P>) -> io::Result<FireResult> {
        if !self.force_tolerance.is_finite()
            || self.force_tolerance <= 0.0
            || !self.initial_dt.is_finite()
            || self.initial_dt <= 0.0
            || !self.max_dt.is_finite()
            || self.max_dt < self.initial_dt
            || !network.box_length.is_finite()
            || network.box_length <= 0.0
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "FIRE requiere tolerancia, pasos de tiempo y caja positivos y finitos; max_dt >= initial_dt",
            ));
        }

        // Las velocidades son auxiliares del minimizador, no estado físico de la red.
        let mut velocities = vec![[0.0; 2]; network.particles.len()];
        let mut dt = self.initial_dt;
        let mut alpha = 0.1;
        let mut positive_steps = 0;
        network.update_forces()?;

        for step in 0..=self.max_steps {
            let mut max_force: f64 = 0.0;
            for particle in &network.particles {
                let force = particle.force[0].hypot(particle.force[1]);
                if !force.is_finite() {
                    return Err(numerical_error());
                }
                max_force = max_force.max(force);
            }
            if max_force <= self.force_tolerance || step == self.max_steps {
                let energy = network.potential_energy();
                if !energy.is_finite() {
                    return Err(numerical_error());
                }
                return Ok(FireResult {
                    converged: max_force <= self.force_tolerance,
                    steps: step,
                    energy,
                    max_force,
                });
            }

            // Velocity Verlet: medio paso de velocidad, posición, nuevas fuerzas.
            for (particle, velocity) in network.particles.iter_mut().zip(&mut velocities) {
                for (axis, component) in velocity.iter_mut().enumerate() {
                    *component += 0.5 * dt * particle.force[axis];
                    let position = particle.position[axis] + dt * *component;
                    if !position.is_finite() {
                        return Err(numerical_error());
                    }
                    particle.position[axis] = position.rem_euclid(network.box_length);
                }
            }
            network.update_forces()?;

            let mut power = 0.0;
            let mut velocity_squared = 0.0;
            let mut force_squared = 0.0;
            for (particle, velocity) in network.particles.iter().zip(&mut velocities) {
                for (axis, component) in velocity.iter_mut().enumerate() {
                    *component += 0.5 * dt * particle.force[axis];
                    power += *component * particle.force[axis];
                    velocity_squared += *component * *component;
                    force_squared += particle.force[axis] * particle.force[axis];
                }
            }
            if !power.is_finite() || !velocity_squared.is_finite() || !force_squared.is_finite() {
                return Err(numerical_error());
            }

            if power > 0.0 {
                // La mezcla usa las normas globales de los 2N grados de libertad.
                let scale = (velocity_squared / force_squared).sqrt();
                for (particle, velocity) in network.particles.iter().zip(&mut velocities) {
                    for (axis, component) in velocity.iter_mut().enumerate() {
                        *component =
                            (1.0 - alpha) * *component + alpha * scale * particle.force[axis];
                    }
                }
                positive_steps += 1;
                if positive_steps > 5 {
                    dt = (dt * 1.1).min(self.max_dt);
                    alpha *= 0.99;
                }
            } else {
                velocities.fill([0.0; 2]);
                positive_steps = 0;
                dt *= 0.5;
                alpha = 0.1;
                if dt == 0.0 {
                    return Err(numerical_error());
                }
            }
        }
        unreachable!("El límite de pasos se comprueba dentro del bucle")
    }
}

fn numerical_error() -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        "FIRE encontró valores numéricos inválidos; revisa el potencial y reduce el paso de tiempo",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::ElasticPotential;
    use std::io::Cursor;

    fn pair() -> Network<ElasticPotential> {
        let data = "0.0 -1.1\n0.1 0.2 1 -1\n0.3 0.2 0 -1\n";
        Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap()
    }

    #[test]
    fn equilibrium_returns_without_moving_particles() {
        let mut network = pair();
        let positions: Vec<_> = network.particles.iter().map(|p| p.position).collect();
        let result = Fire::default().minimize(&mut network).unwrap();
        assert!(result.converged);
        assert_eq!(result.steps, 0);
        assert_eq!(result.energy, 0.0);
        for (particle, original) in network.particles.iter().zip(positions) {
            assert_eq!(particle.position, original);
        }
    }

    #[test]
    fn stretched_pair_relaxes_to_its_rest_length() {
        let mut network = pair();
        network.particles[1].position[0] = 0.4;
        let initial_energy = network.potential_energy();
        let fire = Fire::default();
        let result = fire.minimize(&mut network).unwrap();
        assert!(result.converged);
        assert!(result.steps > 0);
        assert!(result.max_force <= fire.force_tolerance);
        assert!(result.energy < initial_energy);
        assert!(result.energy < 1e-14);
        assert!((network.contact_length(&network.contacts[0]) - 0.2).abs() < 1e-8);
        assert_eq!(network.box_length, 1.0);
    }

    #[test]
    fn affine_expansion_then_relaxation_preserves_box_and_contacts() {
        let data = "0.0 -1.1\n0.1 0.1 1 2\n0.4 0.1 0 2\n0.2 0.4 0 1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        let reference: Vec<_> = network
            .contacts
            .iter()
            .map(|c| (c.particles, c.initial_length))
            .collect();
        network.affine_expand(1.1).unwrap();
        let initial_energy = network.potential_energy();
        let result = Fire::default().minimize(&mut network).unwrap();
        assert!(result.converged);
        assert!(result.energy < initial_energy);
        assert!(result.energy < 1e-14);
        assert_eq!(network.box_length, 1.1);
        assert_eq!(network.contacts.len(), reference.len());
        for (contact, (particles, length)) in network.contacts.iter().zip(reference) {
            assert_eq!(contact.particles, particles);
            assert_eq!(contact.initial_length, length);
        }
    }

    #[test]
    fn periodic_relaxation_wraps_particles_into_the_box() {
        let data = "0.0 -1.1\n0.98 0.5 1 -1\n0.08 0.5 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        network.particles[1].position[0] = 0.28;
        let result = Fire::default().minimize(&mut network).unwrap();
        assert!(result.converged);
        assert!((network.contact_length(&network.contacts[0]) - 0.1).abs() < 1e-8);
        assert!(network.particles[0].position[0] < 0.2);
        assert!(network.particles.iter().all(|p| {
            p.position
                .iter()
                .all(|&x| (0.0..network.box_length).contains(&x))
        }));
    }

    #[test]
    fn equilibrium_does_not_require_zero_energy_or_zero_contact_forces() {
        // Tres resortes alrededor de la caja: las tensiones se compensan.
        let data =
            "0.0 -1.1\n0.0 0.5 1 2\n0.3333333333333333 0.5 0 2\n0.6666666666666666 0.5 0 1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        network.affine_expand(1.1).unwrap();
        let result = Fire::default().minimize(&mut network).unwrap();
        assert!(result.converged);
        assert_eq!(result.steps, 0);
        assert!(result.energy > 0.0);
        assert!(network.contacts.iter().all(|c| c.contact_force > 0.0));
    }

    #[test]
    fn reports_step_limit_without_claiming_convergence() {
        let mut network = pair();
        network.particles[1].position[0] = 0.4;
        let fire = Fire {
            max_steps: 1,
            ..Fire::default()
        };
        let result = fire.minimize(&mut network).unwrap();
        assert!(!result.converged);
        assert_eq!(result.steps, 1);
        assert!(result.max_force > fire.force_tolerance);
        assert_eq!(result.energy, network.potential_energy());
    }

    #[test]
    fn rejects_invalid_parameters_and_propagates_force_errors() {
        for fire in [
            Fire {
                initial_dt: 0.0,
                ..Fire::default()
            },
            Fire {
                max_dt: 0.001,
                ..Fire::default()
            },
            Fire {
                force_tolerance: f64::NAN,
                ..Fire::default()
            },
        ] {
            let mut network = pair();
            assert!(fire.minimize(&mut network).is_err());
            assert_eq!(network.particles[0].position, [0.1, 0.2]);
        }
        let mut network = pair();
        network.particles[1].position = network.particles[0].position;
        assert!(Fire::default().minimize(&mut network).is_err());
    }
}
