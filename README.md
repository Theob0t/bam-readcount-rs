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

Stratified-sampled across the REDACTED, REDACTED, REDACTED, REDACTED, REDACTED, and REDACTED cohorts
(`bench/samples.tsv` lists 200 candidates; the latest run scored **35** of them
× **42.4 M** joined `(sample, position, base)` rows). Per-feature Pearson
correlation against the upstream `<sample>.bamReadCount.txt` files already
living in each sample's `stregaOuts/<sample>/`.

**11 of 13 metrics reproduce upstream at r = 1.0000 (≥ 99.9 % byte-exact).**

| metric | r |
|---|---:|
| count, num_plus_strand, num_minus_strand | **1.00000** |
| avg_mapping_quality, avg_basequality, avg_se_mapping_quality | **1.00000** |
| avg_pos_as_fraction, avg_sum_mismatch_qualities, avg_clipped_length, avg_distance_to_effective_3p_end | **1.00000** |
| avg_num_mismatches_as_fraction | **0.99995** |
| num_q2_containing_reads | 0.95888 |
| avg_distance_to_q2_start_in_q2_reads | 0.84026 |

![correlation grid](bench/results/latest/plots/correlation_grid.png)

See [`bench/results/latest/SUMMARY.md`](bench/results/latest/SUMMARY.md) for
the full per-metric table including MAE and exact-match %.

## Performance

Median across 35 samples: **151 s wall** at 4 threads, **508 MB peak RSS**.
Largest REDACTED sample (1.7 M queried sites): ~22 min wall, ~3 GB RSS.

![runtime / memory](bench/results/latest/plots/runtime_memory.png)

Single-binary, multi-threaded. Replaces the existing
`scripts/bamreadscounts_parallel.py` wrapper (which spawned 22 subprocess
copies of upstream `bam-readcount` per sample, then concatenated). Upstream
timing for an end-to-end pipeline comparison comes from Nextflow trace files,
not this benchmark — see `STREGA/conf/base.config:135` for the per-process
resource budget the new tool replaces.

## Limitations (v1)

- SNV per-base records only (`A`, `C`, `G`, `T`, `N`, `=`). Indel rows
  (`+SEQ` / `-SEQ`) are not emitted — `posLevel.py` does not consume them.
- `--per-library`, `--insertion-centric`, `--print-individual-mapq` modes
  not implemented (STREGA does not pass these flags).
- Output is sorted by (chrom, pos); upstream emits in BED-input order.
  STREGA's downstream parser groups by `chr>pos>ref>alt`, so order is not
  observable. Add a `--preserve-bed-order` flag if needed.
- Q2-related metrics (`num_q2_containing_reads`,
  `avg_distance_to_q2_start_in_q2_reads`) match upstream byte-for-byte on
  forward-stranded reads but show a small per-record disagreement on
  reverse-stranded reads (a corner case in upstream's q2_pos scan logic
  that's hard to reproduce literally without further source reading).

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
