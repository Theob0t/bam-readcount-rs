# Benchmark summary — `200samples_final`

- samples scored: **35**
- joined rows (per-sample × position × base): **42,380,565**
- median rs wall time: **151.3s** (4 threads)
- median rs peak RSS: **508 MB**

## Per-feature correlation

| metric | n | pearson r | MAE | % exact | % within 0.01 |
|---|---:|---:|---:|---:|---:|
| count | 42,380,565 | 1.00000 | 0 | 100.00% | 100.00% |
| avg_mapping_quality | 42,380,565 | 1.00000 | 2.069e-06 | 99.98% | 100.00% |
| avg_basequality | 42,380,565 | 1.00000 | 7.352e-06 | 99.93% | 99.99% |
| avg_se_mapping_quality | 42,380,565 | 1.00000 | 4.766e-08 | 100.00% | 100.00% |
| num_plus_strand | 42,380,565 | 1.00000 | 0 | 100.00% | 100.00% |
| num_minus_strand | 42,380,565 | 1.00000 | 0 | 100.00% | 100.00% |
| avg_pos_as_fraction | 42,380,565 | 1.00000 | 6.706e-06 | 99.93% | 99.94% |
| avg_num_mismatches_as_fraction | 42,380,565 | 0.99995 | 7.386e-06 | 99.93% | 99.99% |
| avg_sum_mismatch_qualities | 42,380,565 | 1.00000 | 7.104e-06 | 99.93% | 99.98% |
| num_q2_containing_reads | 42,380,565 | 0.95888 | 63.31 | 45.89% | 45.89% |
| avg_distance_to_q2_start_in_q2_reads | 42,380,565 | 0.84026 | 0.06879 | 52.72% | 53.62% |
| avg_clipped_length | 42,380,565 | 1.00000 | 7.204e-06 | 99.93% | 99.96% |
| avg_distance_to_effective_3p_end | 42,380,565 | 1.00000 | 2.054e-05 | 99.79% | 99.83% |
