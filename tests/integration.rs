//! Integration tests at the CLI public-interface boundary, per the `tdd`
//! skill: each fixture under `tests/data/fixtures/<name>/` provides a
//! `region.bam` + `region.bed` + `expected.txt` (captured once with upstream
//! `bam-readcount`). The test runs the in-process binary and asserts
//! byte-equal stdout.
//!
//! Add a new fixture by:
//!   1. Carving a tiny BAM region (samtools view -b sample.bam chr:s-e)
//!   2. Writing the BED for the positions you want to query
//!   3. Capturing upstream output once: bam-readcount -w1 -f ref -l region.bed
//!      region.bam > expected.txt
//!   4. Adding the directory name to the FIXTURES list below.

use std::path::PathBuf;
use std::process::Command;

const FIXTURES: &[&str] = &["T01_tracer"];

fn binary_path() -> PathBuf {
    // cargo test sets CARGO_BIN_EXE_<name>; fall back to target/release.
    if let Some(p) = option_env!("CARGO_BIN_EXE_bam-readcount-rs") {
        return PathBuf::from(p);
    }
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("target/release/bam-readcount-rs");
    p
}

fn fasta_path() -> PathBuf {
    PathBuf::from("/gpfs/commons/resources/GatkBundle/Homo_sapiens_assembly38.fasta")
}

#[test]
fn fixtures_match_upstream_byte_for_byte() {
    let bin = binary_path();
    let fa = fasta_path();
    if !bin.exists() {
        eprintln!("skipping: binary not built ({})", bin.display());
        return;
    }
    if !fa.exists() {
        eprintln!("skipping: reference FASTA not present");
        return;
    }

    let mut failures: Vec<String> = Vec::new();
    let fixtures_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/fixtures");

    for name in FIXTURES {
        let dir = fixtures_dir.join(name);
        let bed = dir.join("region.bed");
        let bam = dir.join("region.bam");
        let expected = dir.join("expected.txt");
        if !bed.exists() || !bam.exists() || !expected.exists() {
            failures.push(format!("{name}: fixture files missing"));
            continue;
        }
        let out = dir.join("got.txt");
        let status = Command::new(&bin)
            .arg("-f").arg(&fa)
            .arg("-l").arg(&bed)
            .arg(&bam)
            .arg("-o").arg(&out)
            .status()
            .expect("spawn binary");
        if !status.success() {
            failures.push(format!("{name}: binary exited with {status}"));
            continue;
        }
        let got = std::fs::read_to_string(&out).unwrap();
        let want = std::fs::read_to_string(&expected).unwrap();
        if got != want {
            failures.push(format!("{name}: output diverged from expected.txt"));
        }
    }
    assert!(failures.is_empty(), "fixture mismatches:\n{}", failures.join("\n"));
}
