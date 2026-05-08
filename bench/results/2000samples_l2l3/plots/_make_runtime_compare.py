"""Per-cohort runtime comparison: bam-readcount Original vs bam-readcount-rs.

Bar chart with paired bars (one per tool) for each cohort.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

RUNDIR = Path('/gpfs/commons/home/tbotella/bam-readcount-rs/bench/results/2000samples_l2l3')
TRACE_ROOT = '/gpfs/commons/groups/landau_lab/SCORPIO/PAPER_cohorts'

mpl.rcParams['figure.dpi'] = 110
mpl.rcParams['savefig.dpi'] = 200
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['font.family'] = 'DejaVu Sans'


def parse_dur(s: str) -> float | None:
    if not s or s == '-':
        return None
    total, found = 0.0, False
    for m in re.finditer(r'(\d+(?:\.\d+)?)\s*(ms|s|m|h)', s):
        v, u = float(m.group(1)), m.group(2)
        total += v * {'ms': 1e-3, 's': 1, 'm': 60, 'h': 3600}[u]
        found = True
    return total if found else None


def load_upstream() -> dict[str, float]:
    out: dict[str, float] = {}
    for tf in glob.glob(f'{TRACE_ROOT}/**/trace*.txt', recursive=True):
        try:
            with open(tf) as fh:
                hdr = fh.readline().rstrip('\n').split('\t')
                if not {'name', 'realtime', 'status'} <= set(hdr):
                    continue
                iN, iS, iR = hdr.index('name'), hdr.index('status'), hdr.index('realtime')
                for line in fh:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) <= max(iN, iS, iR):
                        continue
                    if 'BAM_READCOUNTS' not in parts[iN] or parts[iS] not in ('COMPLETED', 'CACHED'):
                        continue
                    m = re.search(r'\(([^)]+)\)', parts[iN])
                    rt = parse_dur(parts[iR])
                    if not m or rt is None:
                        continue
                    sid = m.group(1)
                    if sid not in out or rt > out[sid]:
                        out[sid] = rt
        except OSError:
            pass
    return out


def main() -> None:
    ms = pl.read_csv(RUNDIR / 'per_sample_metrics.tsv', separator='\t').with_columns([
        pl.col('rs_wall_s').cast(pl.Float64),
    ])
    upstream = load_upstream()
    ms = ms.join(
        pl.DataFrame({'sample_id': list(upstream.keys()),
                      'upstream_wall_s': list(upstream.values())}),
        on='sample_id', how='left',
    )
    paired = ms.drop_nulls('upstream_wall_s')

    # Cohort labels are access-controlled — anonymize to cohort_1..cohort_N
    # (deterministic, alphabetical on the underlying name) before plotting.
    real_cohorts = sorted(c for c in paired['cohort'].unique().to_list() if c is not None)
    cohort_anon = {c: f'cohort_{i+1}' for i, c in enumerate(real_cohorts)}
    paired = paired.with_columns(
        pl.col('cohort').replace_strict(cohort_anon).alias('cohort')
    )

    agg = (
        paired.group_by('cohort')
        .agg([
            pl.col('rs_wall_s').median().alias('rs_med'),
            pl.col('upstream_wall_s').median().alias('orig_med'),
            pl.col('rs_wall_s').quantile(0.25).alias('rs_q25'),
            pl.col('rs_wall_s').quantile(0.75).alias('rs_q75'),
            pl.col('upstream_wall_s').quantile(0.25).alias('orig_q25'),
            pl.col('upstream_wall_s').quantile(0.75).alias('orig_q75'),
            pl.len().alias('n'),
        ])
        .sort('orig_med', descending=True)
    )

    cohorts = agg['cohort'].to_list()
    rs_med = agg['rs_med'].to_numpy()
    orig_med = agg['orig_med'].to_numpy()
    rs_err = np.vstack([rs_med - agg['rs_q25'].to_numpy(), agg['rs_q75'].to_numpy() - rs_med])
    orig_err = np.vstack([orig_med - agg['orig_q25'].to_numpy(), agg['orig_q75'].to_numpy() - orig_med])
    ns = agg['n'].to_numpy()

    fig, ax = plt.subplots(figsize=(max(7.5, 1.4 * len(cohorts) + 2), 5.2))
    x = np.arange(len(cohorts))
    w = 0.38
    c_orig, c_rs = '#c44e52', '#3a7ab8'

    b1 = ax.bar(x - w / 2, orig_med, w, yerr=orig_err, capsize=3,
                color=c_orig, edgecolor='black', lw=0.4,
                label='bam-readcount Original', error_kw={'lw': 0.7, 'ecolor': '#333'})
    b2 = ax.bar(x + w / 2, rs_med, w, yerr=rs_err, capsize=3,
                color=c_rs, edgecolor='black', lw=0.4,
                label='bam-readcount-rs', error_kw={'lw': 0.7, 'ecolor': '#333'})

    for xi, (om, rm) in enumerate(zip(orig_med, rs_med)):
        ax.text(xi - w / 2, om, f' {om:.0f}s', ha='center', va='bottom', fontsize=8, rotation=0)
        ax.text(xi + w / 2, rm, f' {rm:.1f}s', ha='center', va='bottom', fontsize=8, rotation=0)
        sp = om / rm if rm > 0 else float('nan')
        ax.text(xi, max(om, rm) * 1.55, f'{sp:.0f}×', ha='center', va='bottom',
                fontsize=9, color='#0a8a3a', fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([f'{c}\n(n={n})' for c, n in zip(cohorts, ns)], fontsize=9)
    ax.set_ylabel('wall time (s, log scale) — median, IQR error bars')
    ax.set_yscale('log')
    ax.set_ylim(top=ax.get_ylim()[1] * 3.5)
    ax.set_title(f'Runtime per cohort — bam-readcount Original vs bam-readcount-rs\n'
                 f'{paired.height} paired samples · overall median speedup '
                 f'{float(np.median(paired["upstream_wall_s"].to_numpy() / paired["rs_wall_s"].to_numpy())):.1f}×')
    ax.legend(frameon=False, loc='upper right')
    ax.grid(axis='y', alpha=0.25, ls='--', lw=0.5)
    ax.set_axisbelow(True)

    out = RUNDIR / 'plots' / 'runtime_compare.png'
    fig.savefig(out)
    print(f'wrote {out}')
    for c, om, rm, n in zip(cohorts, orig_med, rs_med, ns):
        print(f'  {c:<10} n={n:>4}   orig={om:>7.1f}s   rs={rm:>5.1f}s   speedup={om/rm:>5.1f}×')


if __name__ == '__main__':
    main()
