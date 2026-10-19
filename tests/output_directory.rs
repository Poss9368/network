use std::{
    fs,
    path::Path,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

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
    let summary_name = format!("{name}_resumen.dat");
    let output = project.join("output");
    let snapshot = output.join(&snapshot_name);
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
