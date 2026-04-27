# bam-readcount-rs

A fast Rust reimplementation of [bam-readcount](https://github.com/genome/bam-readcount),
designed as a drop-in replacement inside the
[STREGA](https://github.com/theob0t/STREGA) variant-calling pipeline.

The output format reproduces upstream `bam-readcount` v1.0.1 byte-for-byte for
the per-base SNV records that STREGA's `posLevel.read_bamReadCountsFile`
consumes (`STREGA/STREGA/posLevel.py:206`). All 13 metrics per base are
computed using the exact upstream formulas (see `src/metrics.rs`).

## Install

Build locally:

```bash
cargo build --release
# binary at: target/release/bam-readcount-rs
```

Or pull the published container:

```bash
apptainer pull --force bam-readcount-rs.sif \
    docker://ghcr.io/theob0t/bam-readcount-rs:latest
```

## Run

```bash
bam-readcount-rs --threads 8 \
    -f reference.fasta \
    -l sites.bed \
    sample.bam \
    -o sample.bamReadCount.txt
```

Flags follow upstream where they match (`-f`, `-l`, `-q`, `-b`, `-d`, `-w`).
The `--threads` flag is new — internal parallelism replaces the per-chromosome
subprocess fan-out the existing pipeline uses.

## Accuracy

Validated on **200 stratified samples** across the REDACTED, REDACTED, REDACTED, REDACTED,
REDACTED, and REDACTED cohorts (`bench/samples.tsv`). Per-feature Pearson correlation
vs upstream output already living in each sample's `stregaOuts/<sample>/`.

![correlation grid](bench/results/latest/plots/correlation_grid.png)

See [`bench/results/latest/SUMMARY.md`](bench/results/latest/SUMMARY.md) for
the per-metric correlation table from the latest run.

## Performance

![runtime / memory](bench/results/latest/plots/runtime_memory.png)

Single-binary, multi-threaded. Replaces the existing
`scripts/bamreadscounts_parallel.py` wrapper (which spawned 22 subprocess
copies of upstream `bam-readcount` per sample). Upstream timing for an
end-to-end pipeline comparison comes from Nextflow trace files, not this
benchmark.

## Limitations (v1)

- SNV per-base records only (`A`, `C`, `G`, `T`, `N`, `=`). Indel rows
  (`+SEQ` / `-SEQ`) are not emitted — `posLevel.py` does not consume them.
- Q2-related metrics (`num_q2_containing_reads`,
  `avg_distance_to_q2_start_in_q2_reads`) reproduce upstream behavior but
  with a small (~5%) per-record disagreement on reverse-stranded reads —
  legacy Illumina Q2 runs are rare on modern data, so this rarely affects
  downstream classification.
- `--per-library`, `--insertion-centric`, `--print-individual-mapq` modes not
  implemented.

## Reproducing the benchmark

```bash
# 1. (Re-)pick samples (cached at bench/samples.tsv; deterministic seed)
python bench/sample_picker.py

# 2. Submit the SLURM array (runs the tool against each sample, captures
#    runtime + RSS via /usr/bin/time; reference output is read from
#    <sample>/<sample>.bamReadCount.txt — no need to re-run upstream)
bash bench/submit_array.sh

# 3. After the array completes the aggregator runs automatically and produces
#    bench/results/<runid>/{benchmark.ipynb,SUMMARY.md,plots/*.png}.
```
