use crate::potential::Potential;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// La posición inicial se lee en una caja de 1 × 1; luego usa la caja actual.
/// El índice dentro de `Network::particles` es el identificador de la partícula.
#[derive(Debug, Clone)]
pub struct Particle {
    pub position: [f64; 2],
    pub force: [f64; 2],
}

/// Un enlace no dirigido, almacenado una sola vez (primer índice < segundo índice).
#[derive(Debug, Clone)]
pub struct Contact {
    pub particles: [usize; 2],
    /// Distancia inicial con la convención de imagen mínima periódica.
    /// Su interpretación como longitud de reposo dependerá del modelo de fuerzas.
    pub initial_length: f64,
    /// Fuerza escalar con signo: positiva en tensión, negativa en compresión.
    pub contact_force: f64,
}

#[derive(Debug)]
pub struct Network<P: Potential> {
    pub particles: Vec<Particle>,
    pub contacts: Vec<Contact>,
    pub box_length: f64,
    pub potential: P,
}

impl<P: Potential> Network<P> {
    /// Guarda la geometría actual con L en la primera columna de la cabecera.
    pub fn to_file(&self, path: impl AsRef<Path>) -> io::Result<()> {
        let mut writer = BufWriter::new(File::create(path)?);
        self.write_to(&mut writer)?;
        writer.flush()
    }

    /// Mismo formato de vecinos que el input; las posiciones usan la caja actual.
    /// Es una salida para visualizar, no guarda las longitudes de reposo ni el potencial.
    pub fn write_to(&self, mut writer: impl Write) -> io::Result<()> {
        let mut neighbors = vec![Vec::new(); self.particles.len()];
        for contact in &self.contacts {
            let [i, j] = contact.particles;
            if i >= neighbors.len() || j >= neighbors.len() || i == j {
                return Err(invalid_data(
                    "No se puede guardar un contacto con índices inválidos",
                ));
            }
            neighbors[i].push(j);
            neighbors[j].push(i);
        }

        write!(writer, "{:.17e}", self.box_length)?;
        for _ in 0..9 {
            write!(writer, "\t-1.1")?;
        }
        writeln!(writer)?;
        for (particle, adjacent) in self.particles.iter().zip(&mut neighbors) {
            adjacent.sort_unstable();
            write!(
                writer,
                "{:.17e}\t{:.17e}",
                particle.position[0], particle.position[1]
            )?;
            for neighbor in adjacent.iter() {
                write!(writer, "\t{neighbor}")?;
            }
            for _ in adjacent.len()..8 {
                write!(writer, "\t-1")?;
            }
            writeln!(writer)?;
        }
        Ok(())
    }

    /// Lee el formato del generador: cabecera, luego x, y e índices de vecinos.
    pub fn from_file(path: impl AsRef<Path>, potential: P) -> io::Result<Self> {
        Self::from_reader(BufReader::new(File::open(path)?), potential)
    }

    pub fn from_reader(reader: impl BufRead, potential: P) -> io::Result<Self> {
        let mut lines = reader.lines();
        // La primera línea contiene constantes, no una partícula.
        lines
            .next()
            .ok_or_else(|| invalid_data("Falta la cabecera"))??;
        let mut particles = Vec::new();
        let mut neighbors: Vec<Vec<usize>> = Vec::new();

        for (i, line) in lines.enumerate() {
            let line = line?;
            let fields: Vec<_> = line.split_whitespace().collect();
            let line_number = i + 2;
            if fields.len() < 2 {
                return Err(invalid_data(format!("Línea {line_number}: faltan x e y")));
            }

            let mut position = [0.0; 2];
            for axis in 0..2 {
                let value: f64 = fields[axis].parse().map_err(|_| {
                    invalid_data(format!("Línea {line_number}: coordenada inválida"))
                })?;
                if !value.is_finite() || !(0.0..1.0).contains(&value) {
                    return Err(invalid_data(format!(
                        "Línea {line_number}: las coordenadas deben estar en [0, 1)"
                    )));
                }
                position[axis] = value;
            }

            let mut adjacent = Vec::new();
            for field in &fields[2..] {
                if *field == "-1" {
                    continue;
                }
                let neighbor: usize = field.parse().map_err(|_| {
                    invalid_data(format!("Línea {line_number}: índice de vecino inválido"))
                })?;
                if neighbor == i || adjacent.contains(&neighbor) {
                    return Err(invalid_data(format!(
                        "Línea {line_number}: enlace consigo mismo o vecino duplicado"
                    )));
                }
                adjacent.push(neighbor);
            }
            particles.push(Particle {
                position,
                force: [0.0; 2],
            });
            neighbors.push(adjacent);
        }

        if particles.is_empty() {
            return Err(invalid_data("La red no contiene partículas"));
        }
        let mut network = Self {
            particles,
            contacts: Vec::new(),
            box_length: 1.0,
            potential,
        };

        for (i, adjacent) in neighbors.iter().enumerate() {
            for &j in adjacent {
                if j >= network.particles.len() {
                    return Err(invalid_data(format!(
                        "Vecino {j} del nodo {i} fuera de rango"
                    )));
                }
                if !neighbors[j].contains(&i) {
                    return Err(invalid_data(format!("El enlace {i}–{j} no es recíproco")));
                }
                if i < j {
                    let mut contact = Contact {
                        particles: [i, j],
                        initial_length: 0.0,
                        contact_force: 0.0,
                    };
                    contact.initial_length = network.contact_length(&contact);
                    network.contacts.push(contact);
                }
            }
        }
        Ok(network)
    }

    /// Vector de i a la imagen periódica más cercana de j.
    pub fn contact_displacement(&self, contact: &Contact) -> [f64; 2] {
        let [i, j] = contact.particles;
        let [x1, y1] = self.particles[i].position;
        let [x2, y2] = self.particles[j].position;
        let mut dx = x2 - x1;
        let mut dy = y2 - y1;
        dx -= self.box_length * (dx / self.box_length).round();
        dy -= self.box_length * (dy / self.box_length).round();
        [dx, dy]
    }

    pub fn contact_length(&self, contact: &Contact) -> f64 {
        let [dx, dy] = self.contact_displacement(contact);
        dx.hypot(dy)
    }

    pub fn contacts_lengths(&self) -> Vec<f64> {
        self.contacts
            .iter()
            .map(|c| self.contact_length(c))
            .collect()
    }

    /// Energía de la geometría actual, contando cada contacto una sola vez.
    pub fn potential_energy(&self) -> f64 {
        self.contacts
            .iter()
            .map(|contact| {
                self.potential
                    .energy(self.contact_length(contact), contact.initial_length)
            })
            .sum()
    }

    /// Presiones normales [Pxx, Pyy] del virial estático, con área A = L².
    /// Positivas en compresión y negativas en tensión. Usa la geometría y el
    /// potencial actuales, sin depender de fuerzas almacenadas ni de FIRE.
    pub fn normal_pressures(&self) -> io::Result<[f64; 2]> {
        let area = self.box_length * self.box_length;
        if !self.box_length.is_finite()
            || self.box_length <= 0.0
            || !area.is_finite()
            || area == 0.0
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "La presión requiere una caja de lado y área positivos y finitos",
            ));
        }

        let mut pressures = [0.0; 2];
        for contact in &self.contacts {
            let displacement = self.contact_displacement(contact);
            let length = displacement[0].hypot(displacement[1]);
            if !length.is_finite() || length == 0.0 {
                return Err(invalid_data(
                    "La presión requiere contactos de longitud finita y no nula",
                ));
            }
            let force = self.potential.force(length, contact.initial_length);
            if !force.is_finite() {
                return Err(invalid_data("El potencial produjo una fuerza no finita"));
            }
            // force = dU/dl y displacement apunta de i a j: el virial lleva signo menos.
            // Cada contacto se cuenta una sola vez, sin un factor 1/2 adicional.
            for (axis, delta) in displacement.into_iter().enumerate() {
                pressures[axis] -= delta * (force * (delta / length));
            }
        }
        for pressure in &mut pressures {
            *pressure /= area;
            if !pressure.is_finite() {
                return Err(invalid_data("La presión calculada no es finita"));
            }
        }
        Ok(pressures)
    }

    /// Promedio de las presiones normales: P = (Pxx + Pyy) / 2.
    pub fn pressure(&self) -> io::Result<f64> {
        let [px, py] = self.normal_pressures()?;
        Ok(0.5 * px + 0.5 * py)
    }

    /// Fuerzas normales equivalentes por lado: [Pxx * L, Pyy * L].
    /// Son promedios espaciales para la caja periódica, no fuerzas de paredes reales.
    /// El signo positivo es hacia afuera; el negativo, hacia adentro.
    /// Los lados opuestos tienen vectores opuestos, no se suman aquí.
    pub fn wall_forces(&self) -> io::Result<[f64; 2]> {
        let [px, py] = self.normal_pressures()?;
        let forces = [px * self.box_length, py * self.box_length];
        if forces.iter().any(|force| !force.is_finite()) {
            return Err(invalid_data(
                "La fuerza equivalente sobre los lados no es finita",
            ));
        }
        Ok(forces)
    }

    /// Promedio con signo de las fuerzas normales por lado: P * L.
    pub fn mean_wall_force(&self) -> io::Result<f64> {
        let [fx, fy] = self.wall_forces()?;
        Ok(0.5 * fx + 0.5 * fy)
    }

    /// Recalcula las fuerzas desde cero usando la geometría actual.
    /// Si devuelve un error, las fuerzas no deben usarse hasta una actualización válida.
    pub fn update_forces(&mut self) -> io::Result<()> {
        for particle in &mut self.particles {
            particle.force = [0.0; 2];
        }
        for contact in &mut self.contacts {
            contact.contact_force = 0.0;
        }

        for index in 0..self.contacts.len() {
            let contact = &self.contacts[index];
            let [i, j] = contact.particles;
            let displacement = self.contact_displacement(contact);
            let length = displacement[0].hypot(displacement[1]);
            if !length.is_finite() || length == 0.0 {
                return Err(invalid_data(format!(
                    "Contacto {i}–{j}: la longitud debe ser finita y mayor que cero"
                )));
            }
            let force = self.potential.force(length, contact.initial_length);
            if !force.is_finite() {
                return Err(invalid_data(format!(
                    "Contacto {i}–{j}: el potencial produjo una fuerza no finita"
                )));
            }
            self.contacts[index].contact_force = force;
            for (axis, delta) in displacement.into_iter().enumerate() {
                let component = force * (delta / length);
                self.particles[i].force[axis] += component;
                self.particles[j].force[axis] -= component;
            }
        }
        Ok(())
    }

    /// Modifica solo la caja. Para escalar también las posiciones usa `affine_expand`.
    pub fn set_box_length(&mut self, box_length: f64) -> io::Result<()> {
        if !box_length.is_finite() || box_length <= 0.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "El lado de la caja debe ser finito y positivo",
            ));
        }
        self.box_length = box_length;
        Ok(())
    }

    /// Escala la caja y las posiciones, conservando las longitudes de referencia.
    pub fn affine_expand(&mut self, factor: f64) -> io::Result<()> {
        if !factor.is_finite() || factor <= 0.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "El factor debe ser finito y positivo",
            ));
        }
        if self.particles.iter().any(|particle| {
            particle
                .position
                .iter()
                .any(|value| !(value * factor).is_finite())
        }) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "La expansión produce posiciones no finitas",
            ));
        }
        self.set_box_length(self.box_length * factor)?;
        for particle in &mut self.particles {
            particle.position[0] *= factor;
            particle.position[1] *= factor;
        }
        Ok(())
    }

    pub fn mean_coordination(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        2.0 * self.contacts.len() as f64 / self.particles.len() as f64
    }
}

fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::ElasticPotential;
    use std::io::Cursor;

    #[test]
    fn writes_current_box_positions_and_reciprocal_neighbors() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        network.affine_expand(1.5).unwrap();
        let mut output = Vec::new();
        network.write_to(&mut output).unwrap();
        let output = String::from_utf8(output).unwrap();
        let rows: Vec<Vec<_>> = output
            .lines()
            .map(|line| line.split_whitespace().collect())
            .collect();
        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|row| row.len() == 10));
        assert_eq!(rows[0][0].parse::<f64>().unwrap(), 1.5);
        assert!(rows[0][1..].iter().all(|&value| value == "-1.1"));
        for (i, row) in rows[1..].iter().enumerate() {
            assert_eq!(
                row[0].parse::<f64>().unwrap(),
                network.particles[i].position[0]
            );
            assert_eq!(
                row[1].parse::<f64>().unwrap(),
                network.particles[i].position[1]
            );
            assert_eq!(row[2].parse::<usize>().unwrap(), 1 - i);
            assert!(row[3..].iter().all(|&value| value == "-1"));
        }
    }

    #[test]
    fn normal_pressures_distinguish_axes_and_wall_forces_use_current_length() {
        let data = "0.0 -1.1\n0.1 0.1 1 2\n0.3 0.1 0 -1\n0.1 0.4 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        network.affine_expand(1.5).unwrap();
        // Los resortes horizontal y vertical aportan f*l = 0.06 y 0.135.
        let [px, py] = network.normal_pressures().unwrap();
        assert!((px + 0.06 / 2.25).abs() < 1e-12);
        assert!((py + 0.135 / 2.25).abs() < 1e-12);
        assert!((network.pressure().unwrap() - (px + py) / 2.0).abs() < 1e-12);
        let [fx, fy] = network.wall_forces().unwrap();
        assert!((fx + 0.04).abs() < 1e-12);
        assert!((fy + 0.09).abs() < 1e-12);
        assert!((network.mean_wall_force().unwrap() + 0.065).abs() < 1e-12);
        assert!(
            (network.mean_wall_force().unwrap() - network.pressure().unwrap() * network.box_length)
                .abs()
                < 1e-12
        );
        // No se requiere update_forces: la presión usa el estado geométrico actual.
        assert!(
            network
                .contacts
                .iter()
                .all(|contact| contact.contact_force == 0.0)
        );
    }

    #[test]
    fn periodic_pressure_is_positive_in_compression_and_negative_in_tension() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        for (factor, expected_px) in [(1.0, 0.0), (0.5, 0.02), (1.5, -1.0 / 150.0)] {
            let mut network =
                Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
            network.affine_expand(factor).unwrap();
            let [px, py] = network.normal_pressures().unwrap();
            assert!((px - expected_px).abs() < 1e-12);
            assert_eq!(py, 0.0);
            assert!((network.pressure().unwrap() - expected_px / 2.0).abs() < 1e-12);
            // Cambiar el origen no debe cambiar el virial periódico.
            for particle in &mut network.particles {
                particle.position[0] = (particle.position[0] + 0.3 * factor).rem_euclid(factor);
            }
            assert!((network.normal_pressures().unwrap()[0] - expected_px).abs() < 1e-12);
        }
    }

    #[test]
    fn pressure_matches_negative_energy_derivative_with_respect_to_area() {
        let data = "0.0 -1.1\n0.95 0.5 1 2\n0.05 0.5 0 -1\n0.95 0.7 0 -1\n";
        let expanded = |factor| {
            let mut network =
                Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
            network.affine_expand(factor).unwrap();
            network
        };
        let h = 1e-6;
        let center = expanded(1.1);
        let plus = expanded(1.1 + h);
        let minus = expanded(1.1 - h);
        let area_difference = plus.box_length.powi(2) - minus.box_length.powi(2);
        let expected = -(plus.potential_energy() - minus.potential_energy()) / area_difference;
        assert!((center.pressure().unwrap() - expected).abs() < 1e-10);
    }

    #[test]
    fn balanced_particle_forces_can_coexist_with_nonzero_pressure() {
        let data =
            "0.0 -1.1\n0.0 0.5 1 2\n0.3333333333333333 0.5 0 2\n0.6666666666666666 0.5 0 1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        network.affine_expand(1.1).unwrap();
        network.update_forces().unwrap();
        assert!(
            network
                .particles
                .iter()
                .all(|p| p.force[0].hypot(p.force[1]) < 1e-12)
        );
        assert!(network.pressure().unwrap() < 0.0);
        assert!(network.mean_wall_force().unwrap() < 0.0);
    }

    #[test]
    fn pressure_rejects_invalid_box_contacts_and_potential() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            network.box_length = length;
            assert!(network.normal_pressures().is_err());
            assert!(network.pressure().is_err());
            assert!(network.wall_forces().is_err());
            assert!(network.mean_wall_force().is_err());
        }
        network.box_length = 1.0;
        network.potential.stiffness = f64::NAN;
        assert!(network.pressure().is_err());
        network.potential.stiffness = 2.0;
        network.particles[1].position = network.particles[0].position;
        assert!(network.pressure().is_err());
    }

    #[test]
    fn periodic_forces_have_correct_sign_and_do_not_accumulate_between_updates() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        for (factor, expected_force) in [(1.0, 0.0), (1.5, 0.1), (0.5, -0.1)] {
            let mut network =
                Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
            network.affine_expand(factor).unwrap();
            network.particles[0].force = [100.0, -100.0];
            network.contacts[0].contact_force = 100.0;

            // Repetir la actualización debe reemplazar el resultado, no duplicarlo.
            for _ in 0..2 {
                network.update_forces().unwrap();
                assert!((network.contacts[0].contact_force - expected_force).abs() < 1e-12);
                assert!((network.particles[0].force[0] - expected_force).abs() < 1e-12);
                assert!((network.particles[1].force[0] + expected_force).abs() < 1e-12);
                assert_eq!(network.particles[0].force[1], 0.0);
                assert_eq!(network.particles[1].force[1], 0.0);
            }
        }
    }

    #[test]
    fn sums_multiple_contacts_and_conserves_total_force() {
        let data = "0.0 -1.1\n0.1 0.1 1 2\n0.2 0.1 0 -1\n0.1 0.2 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        network.affine_expand(2.0).unwrap();
        network.update_forces().unwrap();
        for (particle, expected) in
            network
                .particles
                .iter()
                .zip([[0.2, 0.2], [-0.2, 0.0], [0.0, -0.2]])
        {
            for (force, expected_component) in particle.force.into_iter().zip(expected) {
                assert!((force - expected_component).abs() < 1e-12);
            }
        }
        for axis in 0..2 {
            let total: f64 = network.particles.iter().map(|p| p.force[axis]).sum();
            assert!(total.abs() < 1e-12);
        }
    }

    #[test]
    fn particle_forces_match_negative_energy_gradient() {
        let data = "0.0 -1.1\n0.1 0.1 1 2\n0.2 0.1 0 -1\n0.1 0.2 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(3.0)).unwrap();
        network.particles[0].position = [0.12, 0.14];
        network.update_forces().unwrap();
        let energy = |network: &Network<ElasticPotential>| -> f64 {
            network
                .contacts
                .iter()
                .map(|contact| {
                    network
                        .potential
                        .energy(network.contact_length(contact), contact.initial_length)
                })
                .sum()
        };
        let h = 1e-6;
        for i in 0..network.particles.len() {
            for axis in 0..2 {
                let original = network.particles[i].position[axis];
                network.particles[i].position[axis] = original + h;
                let plus = energy(&network);
                network.particles[i].position[axis] = original - h;
                let minus = energy(&network);
                network.particles[i].position[axis] = original;
                let expected_force = -(plus - minus) / (2.0 * h);
                assert!((network.particles[i].force[axis] - expected_force).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn rejects_zero_length_and_non_finite_forces() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        let original = network.particles[1].position;
        network.particles[1].position = network.particles[0].position;
        assert!(network.update_forces().is_err());
        assert!(network.particles.iter().all(|p| p.force == [0.0; 2]));

        network.particles[1].position = original;
        network.potential.stiffness = f64::NAN;
        assert!(network.update_forces().is_err());
    }

    #[test]
    fn loads_real_network_and_updates_forces_after_expansion() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("networks/network_4096_4.300_0000.dat");
        let mut network = Network::from_file(path, ElasticPotential::new(1.0)).unwrap();
        assert_eq!(network.particles.len(), 4096);
        assert_eq!(network.contacts.len(), 8806);
        network.update_forces().unwrap();
        assert!(network.particles.iter().all(|p| p.force == [0.0; 2]));

        network.affine_expand(1.01).unwrap();
        network.update_forces().unwrap();
        for contact in &network.contacts {
            assert!((contact.contact_force - 0.01 * contact.initial_length).abs() < 1e-12);
        }
        for axis in 0..2 {
            let total: f64 = network.particles.iter().map(|p| p.force[axis]).sum();
            assert!(total.abs() < 1e-12);
        }
    }

    #[test]
    fn loads_one_contact_and_uses_periodic_distance() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        let network = Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        assert_eq!(network.particles.len(), 2);
        assert_eq!(network.contacts.len(), 1);
        assert_eq!(network.contacts[0].particles, [0, 1]);
        assert_eq!(network.box_length, 1.0);
        assert_eq!(network.particles[0].position, [0.95, 0.5]);
        assert_eq!(network.particles[1].position, [0.05, 0.5]);
        assert!(
            network
                .particles
                .iter()
                .all(|particle| particle.force == [0.0; 2])
        );
        assert!((network.contacts[0].initial_length - 0.1).abs() < 1e-12);
        assert_eq!(network.contacts[0].contact_force, 0.0);
        assert_eq!(network.mean_coordination(), 1.0);
    }

    #[test]
    fn rejects_invalid_topology() {
        for rows in [
            "0.1 0.2 1\n0.3 0.4 -1\n",  // No recíproco.
            "0.1 0.2 2\n0.3 0.4 0\n",   // Fuera de rango.
            "0.1 0.2 0\n",              // Consigo mismo.
            "0.1 0.2 1 1\n0.3 0.4 0\n", // Duplicado.
            "0.1 0.2 -2\n",             // Índice negativo distinto de -1.
            "0.1 0.2 vecino\n",         // Índice no numérico.
        ] {
            let data = format!("0.0 -1.1\n{rows}");
            assert!(Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).is_err());
        }
    }

    #[test]
    fn affine_expansion_scales_periodic_lengths_and_keeps_reference_lengths() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        let initial_length = network.contacts[0].initial_length;
        network.affine_expand(1.5).unwrap();
        assert_eq!(network.box_length, 1.5);
        assert!((network.particles[0].position[0] - 1.425).abs() < 1e-12);
        assert_eq!(network.particles[0].position[1], 0.75);
        assert!((network.contacts_lengths()[0] - 1.5 * initial_length).abs() < 1e-12);
        assert_eq!(network.contacts[0].initial_length, initial_length);
        assert_eq!(network.contacts[0].particles, [0, 1]);
    }

    #[test]
    fn invalid_expansion_does_not_modify_the_network() {
        let data = "0.0 -1.1\n0.95 0.5 1 -1\n0.05 0.5 0 -1\n";
        let mut network =
            Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).unwrap();
        let original = network.particles[0].position;
        for factor in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(network.affine_expand(factor).is_err());
            assert_eq!(network.box_length, 1.0);
            assert_eq!(network.particles[0].position, original);
            assert!(network.set_box_length(factor).is_err());
            assert_eq!(network.box_length, 1.0);
        }
        network.affine_expand(2.0).unwrap();
        assert!(network.affine_expand(f64::MAX).is_err());
        assert_eq!(network.box_length, 2.0);
    }

    #[test]
    fn rejects_invalid_geometry() {
        for data in [
            "",
            "0.0 -1.1\n",
            "0.0 -1.1\nNaN 0.2 -1\n",
            "0.0 -1.1\n1.0 0.2 -1\n",
            "0.0 -1.1\n-0.1 0.2 -1\n",
            "0.0 -1.1\n0.1 inf -1\n",
            "0.0 -1.1\n0.1\n",
        ] {
            assert!(Network::from_reader(Cursor::new(data), ElasticPotential::new(2.0)).is_err());
        }
    }
}
