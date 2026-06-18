#!/usr/bin/env python3
"""
Summarize exome-FG LD from a single merged ld.tsv.gz file (all chroms).

Usage:
    summarize_ld.py [OPTIONS] FILE

Chrom is extracted from the FG_SNP column; no per-chrom files needed.

Outputs (in --outdir, prefixed with {prefix}_r{min_r2}_):
    ld_stats.tsv           per-chromosome summary table
    fig1_variants.png      unique variants per chrom, stacked coding/non-coding
    fig2_pairs.png         pair breakdown per chrom (both/fg-only/ex-only/neither)
    fig3_r2_dist.png       R2 distribution by coding category (violin)
    fig4_coding_frac.png   coding fraction per chrom (line)

NOTE: handles \r\n line endings in input files (transitional, to be removed).
"""

import argparse
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── helpers ──────────────────────────────────────────────────────────────────

def chrom_sort_key(c):
    try:
        return (0, int(c))
    except ValueError:
        return (1, c)



def read_ld(path, nrows=None):
    print(f"Reading {path} ...", file=sys.stderr)
    chunks = []
    pbar = tqdm(unit=' rows', unit_scale=True, desc='Reading', file=sys.stderr)
    for chunk in pd.read_csv(path, sep='\t', chunksize=200_000, nrows=nrows):
        for col in ('is_fg_coding', 'is_ex_coding'):
            if col in chunk.columns and chunk[col].dtype == object:
                chunk[col] = chunk[col].map({'True': True, 'False': False}).astype('bool')
        chunk['R2'] = chunk['R2'].astype('float32')
        chunk['chrom'] = [s.split('_')[0].removeprefix('chr') for s in chunk['FG_SNP']]
        chunks.append(chunk)
        pbar.update(len(chunk))
    pbar.close()
    df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    for col in ('FG_SNP', 'EXOME_SNP'):
        if col in df.columns:
            df[col] = df[col].astype('category')
    return df


# ── stats ─────────────────────────────────────────────────────────────────────

def compute_stats(df, chrom):
    fg = df[['FG_SNP',    'is_fg_coding']].drop_duplicates('FG_SNP')
    ex = df[['EXOME_SNP', 'is_ex_coding']].drop_duplicates('EXOME_SNP')

    n_fg  = len(fg)
    n_ex  = len(ex)
    fg_cod = int(fg['is_fg_coding'].sum())
    ex_cod = int(ex['is_ex_coding'].sum())

    both    = int(( df['is_fg_coding'] &  df['is_ex_coding']).sum())
    fg_only = int(( df['is_fg_coding'] & ~df['is_ex_coding']).sum())
    ex_only = int((~df['is_fg_coding'] &  df['is_ex_coding']).sum())
    neither = int((~df['is_fg_coding'] & ~df['is_ex_coding']).sum())

    return {
        'chrom':                 chrom,
        'n_pairs':               len(df),
        'n_fg_variants':         n_fg,
        'fg_coding':             fg_cod,
        'fg_coding_pct':         round(100 * fg_cod / n_fg, 2) if n_fg else 0,
        'n_ex_variants':         n_ex,
        'ex_coding':             ex_cod,
        'ex_coding_pct':         round(100 * ex_cod / n_ex, 2) if n_ex else 0,
        'n_pairs_both_coding':   both,
        'n_pairs_fg_only':       fg_only,
        'n_pairs_ex_only':       ex_only,
        'n_pairs_neither':       neither,
    }


# ── figures ───────────────────────────────────────────────────────────────────

def fig1_variants(stats, outdir, pfx='', r2_label=''):
    """Stacked bar: unique coding/non-coding variants per chrom, FG and exome."""
    chroms = stats['chrom'].tolist()
    x = np.arange(len(chroms))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, (total_col, cod_col, title) in zip(axes, [
        ('n_fg_variants', 'fg_coding',  'FG variants'),
        ('n_ex_variants', 'ex_coding',  'Exome variants'),
    ]):
        coding    = stats[cod_col].values
        noncoding = stats[total_col].values - coding
        ax.bar(x, noncoding, label='non-coding', color='steelblue')
        ax.bar(x, coding, bottom=noncoding, label='coding', color='tomato')
        ax.set_xticks(x)
        ax.set_xticklabels(chroms, rotation=45, ha='right')
        ax.set_ylabel('Unique variants')
        ax.set_title(title)
        ax.legend()

    fig.suptitle(f'Unique variants in LD (FG ↔ Exome) by chromosome ({r2_label})')
    fig.tight_layout()
    fig.savefig(outdir / f'{pfx}fig1_variants.png', dpi=150)
    plt.close(fig)
    print("Saved fig1_variants.png", file=sys.stderr)


def fig2_pairs(stats, outdir, pfx='', r2_label=''):
    """Stacked bar: pair breakdown (both/fg-only/ex-only/neither) per chrom."""
    chroms = stats['chrom'].tolist()
    x = np.arange(len(chroms))

    cols   = ['n_pairs_neither', 'n_pairs_ex_only', 'n_pairs_fg_only', 'n_pairs_both_coding']
    labels = ['neither', 'exome only', 'FG only', 'both coding']
    colors = ['#aec7e8', '#1f77b4', '#ff7f0e', '#d62728']

    fig, ax = plt.subplots(figsize=(12, 5))
    bottom = np.zeros(len(chroms))
    for col, label, color in zip(cols, labels, colors):
        vals = stats[col].values.astype(float)
        ax.bar(x, vals, bottom=bottom, label=label, color=color)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha='right')
    ax.set_ylabel('Pairs')
    ax.set_title(f'LD pair breakdown by coding status per chromosome ({r2_label})')
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / f'{pfx}fig2_pairs.png', dpi=150)
    plt.close(fig)
    print("Saved fig2_pairs.png", file=sys.stderr)


def fig3_r2_dist(r2_data, outdir, pfx='', r2_label=''):
    """Violin: R2 distribution by coding category, pooled across chroms."""
    keys   = ['both_coding', 'fg_only', 'ex_only', 'neither']
    labels = ['both coding', 'FG only', 'exome only', 'neither']

    present = [(k, l) for k, l in zip(keys, labels) if len(r2_data.get(k, [])) > 0]
    if not present:
        return

    data          = [r2_data[k] for k, _ in present]
    valid_labels  = [l for _, l in present]

    fig, ax = plt.subplots(figsize=(8, 5))
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for pc in parts['bodies']:
        pc.set_alpha(0.7)
    ax.set_xticks(range(1, len(valid_labels) + 1))
    ax.set_xticklabels(valid_labels)
    ax.set_ylabel('R²')
    ax.set_ylim(0, 1)
    ax.set_title(f'R² distribution by coding category — all chromosomes ({r2_label})')
    fig.tight_layout()
    fig.savefig(outdir / f'{pfx}fig3_r2_dist.png', dpi=150)
    plt.close(fig)
    print("Saved fig3_r2_dist.png", file=sys.stderr)


def fig4_coding_frac(stats, outdir, pfx='', r2_label=''):
    """Line plot: coding fraction per chrom for FG and exome."""
    chroms = stats['chrom'].tolist()
    x = np.arange(len(chroms))

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x, stats['fg_coding_pct'], marker='o', label='FG coding %',    color='tomato')
    ax.plot(x, stats['ex_coding_pct'], marker='s', label='Exome coding %', color='steelblue')
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha='right')
    ax.set_ylabel('Coding fraction (%)')
    ax.set_title(f'Coding fraction of variants in LD per chromosome ({r2_label})')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f'{pfx}fig4_coding_frac.png', dpi=150)
    plt.close(fig)
    print("Saved fig4_coding_frac.png", file=sys.stderr)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_file', help='merged ld.tsv.gz file (all chroms concatenated)')
    parser.add_argument('--outdir', default='.', help='Output directory (default: cwd)')
    parser.add_argument('--prefix', default='exome_ld', help='Output file prefix (default: exome_ld)')
    parser.add_argument('--min_r2', type=float, default=0.0,
                        help='Keep only pairs with R2 >= threshold (default: 0, no filter)')
    parser.add_argument('--test', nargs='?', const=1000, default=None, type=int,
                        help='Only read first N lines per file (default N=1000 if flag given without value)')
    parser.add_argument('--r2-sample', type=int, default=50_000,
                        help='Max R2 values per category for violin plot (default: 50000)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pfx = f"{args.prefix}_r{args.min_r2}_" if args.prefix else f"r{args.min_r2}_"
    r2_label = f"R²≥{args.min_r2}"

    filtered_path = outdir / f'{pfx}ld.tsv.gz'
    cache_path    = outdir / f'{pfx}df.pkl'

    if cache_path.exists():
        print(f"Loading cached dataframe from {cache_path} ...", file=sys.stderr)
        df = pd.read_pickle(cache_path)
        mem_gb = df.memory_usage(deep=True).sum() / 1e9
        print(f"Loaded {len(df):,} rows from cache — {mem_gb:.2f} GB in memory", file=sys.stderr)
    else:
        if not filtered_path.exists():
            cmd = (f"zcat {shlex.quote(str(args.input_file))}"
                   f" | awk 'NR==1 || $3 >= {args.min_r2}'"
                   f" | bgzip -@ {os.cpu_count()} > {shlex.quote(str(filtered_path))}")
            print(f"Filtering with awk: {cmd}", file=sys.stderr)
            subprocess.run(cmd, shell=True, check=True)
            print(f"Filtered pairs written to {filtered_path}", file=sys.stderr)
        df = read_ld(filtered_path, nrows=args.test)
        mem_gb = df.memory_usage(deep=True).sum() / 1e9
        print(f"Loaded {len(df):,} rows — {mem_gb:.2f} GB in memory", file=sys.stderr)
        print(f"Caching dataframe to {cache_path} ...", file=sys.stderr)
        df.to_pickle(cache_path)

    groups = {chrom: cdf for chrom, cdf in df.groupby('chrom', sort=False)}
    chroms = sorted(groups.keys(), key=chrom_sort_key)
    print(f"Chroms found: {chroms}", file=sys.stderr)

    stats_rows = []
    r2_data = {k: [] for k in ('both_coding', 'fg_only', 'ex_only', 'neither')}
    rng = np.random.default_rng(seed=42)

    for i, chrom in enumerate(chroms):
        print(f"Processing chrom {chrom} ({i+1}/{len(chroms)}) ...", file=sys.stderr)
        cdf = groups[chrom]
        stats_rows.append(compute_stats(cdf, chrom))

        masks = {
            'both_coding': ( cdf['is_fg_coding'] &  cdf['is_ex_coding']),
            'fg_only':     ( cdf['is_fg_coding'] & ~cdf['is_ex_coding']),
            'ex_only':     (~cdf['is_fg_coding'] &  cdf['is_ex_coding']),
            'neither':     (~cdf['is_fg_coding'] & ~cdf['is_ex_coding']),
        }
        for k, mask in masks.items():
            vals = cdf.loc[mask, 'R2'].values
            if len(vals) > args.r2_sample:
                vals = rng.choice(vals, args.r2_sample, replace=False)
            r2_data[k].extend(vals.tolist())

    stats = pd.DataFrame(stats_rows)

    stats_path = outdir / f'{pfx}ld_stats.tsv'
    stats.to_csv(stats_path, sep='\t', index=False)
    print(f"Stats written to {stats_path}", file=sys.stderr)

    # cap total R2 sample across chroms for violin
    rng = np.random.default_rng(seed=42)
    r2_capped = {}
    for k, vals in r2_data.items():
        arr = np.array(vals)
        cap = args.r2_sample * 4
        if len(arr) > cap:
            arr = rng.choice(arr, cap, replace=False)
        r2_capped[k] = arr

    print("Generating fig1_variants ...", file=sys.stderr)
    fig1_variants(stats, outdir, pfx, r2_label)
    print("Generating fig2_pairs ...", file=sys.stderr)
    fig2_pairs(stats, outdir, pfx, r2_label)
    print("Generating fig3_r2_dist ...", file=sys.stderr)
    fig3_r2_dist(r2_capped, outdir, pfx, r2_label)
    print("Generating fig4_coding_frac ...", file=sys.stderr)
    fig4_coding_frac(stats, outdir, pfx, r2_label)

    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
