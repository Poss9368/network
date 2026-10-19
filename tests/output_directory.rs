use std::{
    fs,
    path::Path,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

#[test]
fn no_snapshots_keeps_the_same_summary_with_or_without_custom_first_length() {
    let project = Path::new(env!("CARGO_MANIFEST_DIR"));
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let name = format!("summary-only-test-{}-{stamp}", std::process::id());
    let input = std::env::temp_dir().join(format!("{name}.dat"));
    fs::write(
        &input,
        "0.0 -1.1\n0.0 0.5 1 2\n0.3333333333333333 0.5 0 2\n0.6666666666666666 0.5 0 1\n",
    )
    .unwrap();
    let output = project.join("output");
    let summary = output.join(format!("{name}_summary.dat"));
    for args in [vec!["1.1", "3"], vec!["1.1", "3", "1.001"]] {
        let saved = Command::new(env!("CARGO_BIN_EXE_network"))
            .arg(&input)
            .args(&args)
            .output()
            .unwrap();
        assert!(
            saved.status.success(),
            "{}",
            String::from_utf8_lossy(&saved.stderr)
        );
        let expected = fs::read_to_string(&summary).unwrap();
        assert_eq!(expected.lines().count(), 4);
        for step in 1..=3 {
            fs::remove_file(
                output
                    .join("steps_details")
                    .join(format!("{name}_step_{step:04}.dat")),
            )
            .unwrap();
        }
        fs::remove_file(&summary).unwrap();

        let summary_only = Command::new(env!("CARGO_BIN_EXE_network"))
            .arg(&input)
            .args(&args)
            .arg("--no-snapshots")
            .output()
            .unwrap();
        assert!(
            summary_only.status.success(),
            "{}",
            String::from_utf8_lossy(&summary_only.stderr)
        );
        assert_eq!(fs::read_to_string(&summary).unwrap(), expected);
        for step in 1..=3 {
            assert!(
                !output
                    .join("steps_details")
                    .join(format!("{name}_step_{step:04}.dat"))
                    .exists()
            );
        }
        let files: Vec<_> = fs::read_dir(&output)
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .filter(|file| file.to_string_lossy().starts_with(&name))
            .collect();
        assert_eq!(files, vec![summary.file_name().unwrap().to_os_string()]);
        fs::remove_file(&summary).unwrap();
    }
    fs::remove_file(input).unwrap();
}

#[test]
fn command_accepts_custom_first_length_and_writes_requested_schedule() {
    let project = Path::new(env!("CARGO_MANIFEST_DIR"));
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let name = format!("custom-first-length-{}-{stamp}", std::process::id());
    let input = std::env::temp_dir().join(format!("{name}.dat"));
    fs::write(
        &input,
        "0.0 -1.1\n0.0 0.5 1 2\n0.3333333333333333 0.5 0 2\n0.6666666666666666 0.5 0 1\n",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_network"))
        .arg(&input)
        .args(["1.1", "3", "1.001"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let output = project.join("output");
    let summary = output.join(format!("{name}_summary.dat"));
    let records = fs::read_to_string(&summary).unwrap();
    let lengths: Vec<f64> = records
        .lines()
        .skip(1)
        .map(|line| line.split_whitespace().next().unwrap().parse().unwrap())
        .collect();
    assert_eq!(lengths.len(), 3);
    for (index, expected) in [1.001, 1.01, 1.1].into_iter().enumerate() {
        assert!((lengths[index] - expected).abs() < 1e-14);
        let snapshot = output
            .join("steps_details")
            .join(format!("{name}_step_{:04}.dat", index + 1));
        let data = fs::read_to_string(&snapshot).unwrap();
        let actual: f64 = data.split_whitespace().next().unwrap().parse().unwrap();
        assert_eq!(actual, lengths[index]);
        fs::remove_file(snapshot).unwrap();
    }
    fs::remove_file(summary).unwrap();
    fs::remove_file(input).unwrap();
}

#[test]
fn running_from_networks_writes_to_project_output() {
    let project = Path::new(env!("CARGO_MANIFEST_DIR"));
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let name = format!("output-location-test-{}-{stamp}", std::process::id());
    let temporary = std::env::temp_dir().join(&name);
    fs::create_dir_all(&temporary).unwrap();
    let input = temporary.join(format!("{name}.dat"));
    fs::write(&input, "0.0 -1.1\n0.1 0.2 1 -1\n0.3 0.2 0 -1\n").unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_network"))
        .current_dir(project.join("networks"))
        .arg(&input)
        .args(["1.01", "1"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let snapshot_name = format!("{name}_step_0001.dat");
    let summary_name = format!("{name}_summary.dat");
    let output = project.join("output");
    let snapshot = output.join("steps_details").join(&snapshot_name);
    assert!(!output.join(&snapshot_name).exists());
    let summary = output.join(&summary_name);
    assert!(snapshot.is_file());
    assert!(summary.is_file());
    assert!(
        !project
            .join("networks/output")
            .join(&snapshot_name)
            .exists()
    );
    assert!(!project.join("networks/output").join(&summary_name).exists());
    let data = fs::read_to_string(&snapshot).unwrap();
    assert_eq!(
        data.split_whitespace()
            .next()
            .unwrap()
            .parse::<f64>()
            .unwrap(),
        1.01
    );
    let records = fs::read_to_string(&summary).unwrap();
    assert_eq!(records.lines().count(), 2);

    fs::remove_file(snapshot).unwrap();
    fs::remove_file(summary).unwrap();
    fs::remove_dir_all(temporary).unwrap();
}
