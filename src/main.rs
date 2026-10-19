pub mod fire;
pub mod network;
pub mod potential;

use fire::Fire;
use network::Network;
use potential::ElasticPotential;
use std::{
    env,
    error::Error,
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::ExitCode,
};

/// Rango log10 del strain positivo ε = L - 1, desde ε_final / 1000.
const STRAIN_DECADES: f64 = 3.0;

fn logarithmic_lengths(final_length: f64, steps: usize) -> Result<Vec<f64>, Box<dyn Error>> {
    if !final_length.is_finite() || final_length <= 1.0 {
        return Err("El largo final debe ser finito y mayor que 1".into());
    }
    if steps == 0 {
        return Err("La cantidad de pasos debe ser al menos 1".into());
    }
    let final_strain = final_length - 1.0;
    let mut lengths = Vec::new();
    let mut previous = 1.0;
    for index in 0..steps {
        let length = if index == steps - 1 {
            final_length
        } else {
            let fraction = index as f64 / (steps - 1) as f64;
            1.0 + final_strain * 10.0_f64.powf(STRAIN_DECADES * (fraction - 1.0))
        };
        if length <= previous {
            return Err("Los largos de los pasos no se pueden distinguir con precisión f64; aumenta el strain final o reduce los pasos".into());
        }
        lengths.push(length);
        previous = length;
    }
    Ok(lengths)
}

fn stretch_network(path: &Path, lengths: &[f64], output: &Path) -> Result<(), Box<dyn Error>> {
    let stem = path
        .file_stem()
        .ok_or("El archivo debe tener un nombre válido")?;
    let mut network = Network::from_file(path, ElasticPotential::new(1.0))?;
    fs::create_dir_all(output)?;
    let mut summary_name = stem.to_os_string();
    summary_name.push("_resumen.dat");
    let summary_path = output.join(summary_name);
    let mut summary = BufWriter::new(File::create(&summary_path)?);
    writeln!(summary, "# box_length\tpressure")?;
    let fire = Fire::default();

    for (index, &length) in lengths.iter().enumerate() {
        let step = index + 1;
        // Cada objetivo es un largo absoluto, no un factor respecto del paso anterior.
        network.affine_expand(length / network.box_length)?;
        network.set_box_length(length)?;
        let result = fire.minimize(&mut network)?;
        if !result.converged {
            return Err(format!(
                "FIRE no convergió en el paso {step}; fuerza máxima = {:.6e}",
                result.max_force
            )
            .into());
        }
        let pressure = network.pressure()?;
        let mut name = stem.to_os_string();
        name.push(format!("_step_{step:04}.dat"));
        network.to_file(output.join(name))?;
        writeln!(summary, "{length:.17e}\t{pressure:.17e}")?;
        summary.flush()?;
        println!(
            "Paso {step}/{}: L = {length:.10e}, P = {pressure:.10e}, FIRE = {} pasos",
            lengths.len(),
            result.steps
        );
    }
    println!("Resumen: {}", summary_path.display());
    Ok(())
}

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os().skip(1);
    let mut path = PathBuf::from(
        args.next()
            .ok_or("Uso: cargo run --release -- <archivo.dat> <largo_final> <pasos>")?,
    );
    let final_length: f64 = args
        .next()
        .ok_or("Falta el largo final, que debe ser mayor que 1")?
        .to_string_lossy()
        .parse()?;
    let steps: usize = args
        .next()
        .ok_or("Falta la cantidad de pasos")?
        .to_string_lossy()
        .parse()?;
    if args.next().is_some() {
        return Err(
            "Se esperan exactamente el nombre de la red, largo final y cantidad de pasos".into(),
        );
    }
    let lengths = logarithmic_lengths(final_length, steps)?;
    let project_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    if !path.is_file() && path.components().count() == 1 {
        path = project_dir.join("networks").join(path);
    }
    stretch_network(&path, &lengths, &project_dir.join("output"))
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("Error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn strain_is_spaced_logarithmically_and_reaches_exact_final_length() {
        let lengths = logarithmic_lengths(1.1, 4).unwrap();
        for (actual, expected) in lengths.iter().zip([1.0001, 1.001, 1.01, 1.1]) {
            assert!((actual - expected).abs() < 1e-14);
        }
        for pair in lengths.windows(2) {
            let log_step = (pair[1] - 1.0).log10() - (pair[0] - 1.0).log10();
            assert!((log_step - 1.0).abs() < 1e-10);
        }
        assert_eq!(lengths.last().copied(), Some(1.1));
        assert_eq!(logarithmic_lengths(1.2, 1).unwrap(), vec![1.2]);
    }

    #[test]
    fn rejects_invalid_or_unrepresentable_schedules() {
        for length in [1.0, 0.5, f64::NAN, f64::INFINITY] {
            assert!(logarithmic_lengths(length, 4).is_err());
        }
        assert!(logarithmic_lengths(1.1, 0).is_err());
        assert!(logarithmic_lengths(1.0 + f64::EPSILON, 4).is_err());
    }

    #[test]
    fn writes_relaxed_steps_and_matching_length_pressure_summary() {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = env::temp_dir().join(format!(
            "network-output-test-{}-{stamp}",
            std::process::id()
        ));
        fs::create_dir_all(&directory).unwrap();
        let input = directory.join("sample.dat");
        let output = directory.join("output");
        // Anillo periódico: fuerzas equilibradas y presión analítica no nula.
        fs::write(
            &input,
            "0.0 -1.1\n0.0 0.5 1 2\n0.3333333333333333 0.5 0 2\n0.6666666666666666 0.5 0 1\n",
        )
        .unwrap();
        let lengths = logarithmic_lengths(1.01, 3).unwrap();
        stretch_network(&input, &lengths, &output).unwrap();
        let summary = fs::read_to_string(output.join("sample_resumen.dat")).unwrap();
        let records: Vec<Vec<f64>> = summary
            .lines()
            .skip(1)
            .map(|line| {
                line.split_whitespace()
                    .map(|value| value.parse().unwrap())
                    .collect()
            })
            .collect();
        assert_eq!(records.len(), lengths.len());
        assert_eq!(fs::read_dir(&output).unwrap().count(), lengths.len() + 1);
        for (index, &length) in lengths.iter().enumerate() {
            assert_eq!(records[index].len(), 2);
            assert_eq!(records[index][0], length);
            let expected_pressure = -(length - 1.0) / (6.0 * length);
            assert!((records[index][1] - expected_pressure).abs() < 1e-12);
            let snapshot =
                fs::read_to_string(output.join(format!("sample_step_{:04}.dat", index + 1)))
                    .unwrap();
            let rows: Vec<Vec<f64>> = snapshot
                .lines()
                .map(|line| {
                    line.split_whitespace()
                        .map(|value| value.parse().unwrap())
                        .collect()
                })
                .collect();
            assert_eq!(rows.len(), 4);
            assert!(rows.iter().all(|row| row.len() == 10));
            assert_eq!(rows[0][0], length);
            for (i, row) in rows[1..].iter().enumerate() {
                assert!((row[0] - i as f64 * length / 3.0).abs() < 1e-12);
                assert!((row[1] - 0.5 * length).abs() < 1e-12);
                let expected: Vec<_> = (0..3).filter(|&j| j != i).map(|j| j as f64).collect();
                assert_eq!(row[2..4], expected);
                assert!(row[4..].iter().all(|&value| value == -1.0));
            }
        }
        fs::remove_dir_all(directory).unwrap();
    }
}
