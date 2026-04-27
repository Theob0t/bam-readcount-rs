"""Stream-parse a bam-readcount text file into a long-form parquet.

Memory-efficient pre-processing for the benchmark notebook: each sample's
ref+rs outputs get parsed once into a parquet file, joined on (chrom, pos,
base), and stored as bench/results/<runid>/joined/<sample>.parquet. The
notebook then scans those parquets via polars (lazy, low RAM).

Usage:
    python parse_brc.py <ref_brc> <rs_brc> <out_parquet>
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

INFO = [
    "count",
    "avg_mapping_quality",
    "avg_basequality",
    "avg_se_mapping_quality",
    "num_plus_strand",
    "num_minus_strand",
    "avg_pos_as_fraction",
    "avg_num_mismatches_as_fraction",
    "avg_sum_mismatch_qualities",
    "num_q2_containing_reads",
    "avg_distance_to_q2_start_in_q2_reads",
    "avg_clipped_length",
    "avg_distance_to_effective_3p_end",
]


def parse_brc(path: Path) -> pl.DataFrame:
    """Stream-parse one bam-readcount text file into a long-form polars frame.
    One row per (chrom, pos, base) where count > 0.
    """
    chroms, poss, bases = [], [], []
    cols = {c: [] for c in INFO}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            for field in parts[4:]:
                vals = field.split(":")
                base = vals[0]
                if base == "=":
                    continue
                if int(vals[1]) == 0:
                    continue
                chroms.append(chrom)
                poss.append(pos)
                bases.append(base)
                for k, v in zip(INFO, vals[1:]):
                    cols[k].append(float(v))
    df = pl.DataFrame(
        {
            "chrom": chroms,
            "pos": poss,
            "base": bases,
            **cols,
        }
    )
    return df


def main():
    ref_path, rs_path, out_path = (
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        Path(sys.argv[3]),
    )
    ref_df = parse_brc(ref_path).rename({c: f"ref_{c}" for c in INFO})
    rs_df = parse_brc(rs_path).rename({c: f"rs_{c}" for c in INFO})
    joined = ref_df.join(rs_df, on=["chrom", "pos", "base"], how="inner")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    joined.write_parquet(out_path, compression="zstd")
    print(
        f"{ref_path.name}: ref={len(ref_df):,} rs={len(rs_df):,} joined={len(joined):,} -> {out_path}"
    )


if __name__ == "__main__":
    main()
