# Benchmark summary — `2000samples`

- samples scored: **1955**
- joined rows (per-sample × position × base): **2,316,674,958**
- median rs wall time: **67.2s** (4 threads)
- median rs peak RSS: **620 MB**

## Per-feature correlation

| metric | n | pearson r | MAE | % exact | % within 0.01 |
|---|---:|---:|---:|---:|---:|
| count | 97,750 | 1.00000 | 0 | 100.0000% | 100.0000% |
| avg_mapping_quality | 97,750 | 1.00000 | 1.535e-06 | 99.9847% | 100.0000% |
| avg_basequality | 97,750 | 1.00000 | 8.696e-06 | 99.9130% | 99.9928% |
| avg_se_mapping_quality | 97,750 | 1.00000 | 7.161e-07 | 99.9928% | 100.0000% |
| num_plus_strand | 97,750 | 1.00000 | 0 | 100.0000% | 100.0000% |
| num_minus_strand | 97,750 | 1.00000 | 0 | 100.0000% | 100.0000% |
| avg_pos_as_fraction | 97,750 | 1.00000 | 1.381e-05 | 99.8619% | 99.8691% |
| avg_num_mismatches_as_fraction | 97,750 | 0.99991 | 1.453e-05 | 99.8547% | 99.9898% |
| avg_sum_mismatch_qualities | 97,750 | 0.99964 | 0.01069 | 99.8793% | 99.9130% |
| num_q2_containing_reads | 97,750 | 0.95988 | 38.42 | 52.9974% | 52.9974% |
| avg_distance_to_q2_start_in_q2_reads | 97,750 | 0.85456 | 0.06305 | 57.8261% | 58.5003% |
| avg_clipped_length | 97,750 | 1.00000 | 1.043e-05 | 99.8957% | 99.9294% |
| avg_distance_to_effective_3p_end | 97,750 | 0.99999 | 5.033e-05 | 99.5018% | 99.5632% |
