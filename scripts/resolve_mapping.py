#!/usr/bin/env python3
"""resolve_mapping.py — build final QRY→REF mapping from combined_summary.tsv.

Groups: MATCHED (ID_CONFIRMED, RESOLVED_BY_ID, RESOLVED_BY_ALIAS,
        INFERRED_BY_ELIMINATION, UNIQUE, CONFLICT_KEPT) |
        DROPPED (CONFLICT_DROPPED, AMBIGUOUS_*) | NO MATCH (MISSING)

Usage:
  python resolve_mapping.py combined_summary.tsv [--aliases FILE] [--seed N] [--out FILE]
  python resolve_mapping.py --test [--seed N]
"""
import argparse, os, pathlib, random, sys
import pandas as pd

# test/ folder sits one level above this script (repo root/test/)
_SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
TEST_DIR     = _SCRIPT_DIR.parent / "test"
TEST_SUMMARY = TEST_DIR / "combined_summary.tsv"
TEST_ALIASES = TEST_DIR / "aliases.tsv"

RESOLVED = {"ID_CONFIRMED","RESOLVED_BY_ID","RESOLVED_BY_ALIAS","UNIQUE"}
CONFLICT_PRIORITY = {"ID_CONFIRMED":0,"RESOLVED_BY_ALIAS":1,"RESOLVED_BY_ID":2,"UNIQUE":3}
DEFAULT_ALIASES = "/mnt/disks/data/samples/exclusions/finngen_R12_duplicate_list.txt"

SUMMARY_GROUPS = [
    ("MATCHED", "samples with a final QRY→REF mapping in the output", [
        ("ID_CONFIRMED",    "single candidate; KING match confirmed by matching IDs"),
        ("RESOLVED_BY_ID",  "twins in ref; query ID matched one candidate"),
        ("RESOLVED_BY_ALIAS","twins in ref; candidates are known aliases of each other"),
        ("UNIQUE",          "single candidate; matched by genetics only"),
        ("CONFLICT_KEPT",   "contested ref ID; kept after priority tiebreak"),
    ]),
    ("DROPPED", "found by KING but excluded from final mapping", [
        ("CONFLICT_DROPPED",     "contested ref ID; lost tiebreak; REF_MAPPED = NA"),
        ("AMBIGUOUS_UNRESOLVED", "multiple ref candidates; no resolution possible"),
    ]),
    ("NO MATCH", "absent from ref or below KING concordance threshold", [
        ("MISSING", "no KING match found"),
    ]),
]

# ── Alias loading ─────────────────────────────────────────────────────────────

def load_aliases(path):
    if not path or not os.path.exists(path):
        return {}
    return _parse_alias_lines(open(path), label=path)


def _parse_alias_lines(lines, label=""):
    groups, n = {}, 0
    for line in lines:
        ids = [x.strip() for x in line.strip().split("\t") if x.strip()]
        if len(ids) < 2:
            continue
        g = frozenset(ids)
        for id_ in ids:
            groups[id_] = g   # maps every ID → its alias frozenset
        n += 1
    if n:
        print(f"Loaded {n} alias groups ({len(groups)} IDs) from {label}  "
              f"e.g. {next(iter(groups))!r}")
    return groups


# ── Alias resolution ──────────────────────────────────────────────────────────

def _alias_resolve(query, candidates, ag):
    """Return (ref_candidate | None, note).

    Alias handling is intentionally QRY-only: only the QUERY is looked up in
    the alias file.  REF candidates are never cross-referenced — they are the
    authoritative identities.  The question asked is purely:
        'which of these REF candidates is listed as a known alias of this query?'
    """
    alias_group = ag.get(query)
    if alias_group is None:
        return None, "query_not_in_alias_file"

    # alias_group is the frozenset of all IDs known to be the same person as query
    m = [c for c in candidates if c in alias_group]  # REF candidates that match

    if len(m) == 1:
        return m[0], ""
    if len(m) == 0:
        return None, "no_candidates_in_alias_group"
    return None, f"{len(m)}_candidates_in_alias_group"


def _collect_ref_ids(df):
    """All IDs that appear as candidates in the DUPLICATES column (= known REF IDs)."""
    ids = set()
    for raw in df["DUPLICATES"].dropna():
        raw = str(raw).strip()
        if raw in ("MISSING", "nan", ""):
            continue
        for c in raw.split(","):
            c = c.strip()
            if c:
                ids.add(c)
    return ids


# ── Core pipeline ─────────────────────────────────────────────────────────────

def _row(dataset, query, ref, cands, status, note=""):
    return dict(DATASET=dataset, QUERY=query, REF_MAPPED=ref,
                CANDIDATES=cands, STATUS=status, ALIAS_NOTE=note)


def initial_categorise(df, ag=None):
    ag = ag or {}
    records = []
    for _, row in df.iterrows():
        dataset, query, raw = row["DATASET"], str(row["QUERY"]).strip(), str(row["DUPLICATES"]).strip()
        if raw in ("MISSING", "nan", ""):
            records.append(_row(dataset, query, "NA", "", "MISSING")); continue
        seen, cands = set(), []
        for c in raw.split(","):
            c = c.strip()
            if c and c not in seen: cands.append(c); seen.add(c)
        if len(cands) == 1:
            ref = cands[0]
            if query == ref:
                records.append(_row(dataset, query, ref, raw, "ID_CONFIRMED"))
            elif ag:
                resolved, _ = _alias_resolve(query, cands, ag)
                records.append(_row(dataset, query, ref, raw,
                                    "RESOLVED_BY_ALIAS" if resolved else "UNIQUE"))
            else:
                records.append(_row(dataset, query, ref, raw, "UNIQUE"))
        else:
            if query in cands:
                records.append(_row(dataset, query, query, raw, "RESOLVED_BY_ID"))
            elif ag:
                resolved, note = _alias_resolve(query, cands, ag)
                records.append(_row(dataset, query, resolved or "AMBIGUOUS", raw,
                                    "RESOLVED_BY_ALIAS" if resolved else "AMBIGUOUS_UNRESOLVED", note))
            else:
                records.append(_row(dataset, query, "AMBIGUOUS", raw, "AMBIGUOUS_UNRESOLVED"))
    return pd.DataFrame(records)



def check_surjectivity(result, rng):
    ref_to_idx = {}
    for idx, row in result[result["STATUS"].isin(RESOLVED)].iterrows():
        ref_to_idx.setdefault(row["REF_MAPPED"], []).append(idx)
    for ref_id, indices in ref_to_idx.items():
        if len(indices) == 1: continue
        def pri(i): return CONFLICT_PRIORITY.get(result.at[i, "STATUS"].split("[")[0], 99)
        best = min(pri(i) for i in indices)
        winner = rng.choice([i for i in indices if pri(i) == best])
        for idx in indices:
            orig = result.at[idx, "STATUS"]
            if idx == winner: result.at[idx, "STATUS"] = f"CONFLICT_KEPT[{orig}]"
            else:             result.at[idx, "STATUS"] = f"CONFLICT_DROPPED[{orig}]"; result.at[idx, "REF_MAPPED"] = "NA"
    return result


# ── Stats table ───────────────────────────────────────────────────────────────

def _pct(n, t): return f"{n/t*100:.1f}%" if t else "n/a"


def build_stats_table(result):
    n, datasets = len(result), sorted(result["DATASET"].unique())
    gc = result["STATUS"].value_counts()
    dsc = {ds: result[result["DATASET"] == ds]["STATUS"].value_counts() for ds in datasets}
    def pfx(p, vc): return sum(v for k, v in vc.items() if k.startswith(p))
    cols = ["SECTION","GROUP","STATUS","TOTAL"] + datasets + ["PCT","NOTES"]
    def blank(): return {c:"" for c in cols}
    tot, brk = [], []
    for gh, gn, statuses in SUMMARY_GROUPS:
        prefixes = [s for s,_ in statuses]
        gcount = sum(pfx(p, gc) for p in prefixes)
        row = dict(SECTION="TOTALS", GROUP=gh, STATUS="", TOTAL=gcount, PCT=_pct(gcount, n), NOTES=gn)
        for ds in datasets: row[ds] = sum(pfx(p, dsc[ds]) for p in prefixes)
        tot.append(row)
        for prefix, desc in statuses:
            cnt = pfx(prefix, gc)
            if cnt == 0: continue
            if prefix == "CONFLICT_KEPT":
                kept = result[result["STATUS"].str.startswith("CONFLICT_KEPT")]
                nc = kept["REF_MAPPED"].nunique()
                desc += f"; {nc} ref IDs contested, avg {gcount/nc:.1f} queries/ref" if nc else ""
            row = dict(SECTION="BREAKDOWN", GROUP=gh, STATUS=prefix, TOTAL=cnt, PCT=_pct(cnt, n), NOTES=desc)
            for ds in datasets: row[ds] = pfx(prefix, dsc[ds])
            brk.append(row)
    grand = dict(SECTION="TOTALS", GROUP="TOTAL", STATUS="", TOTAL=n, PCT="100.0%", NOTES="")
    for ds in datasets: grand[ds] = (result["DATASET"] == ds).sum()
    tot.append(grand)
    sep = blank(); sep["SECTION"] = "---"
    return pd.DataFrame(tot + [sep] + brk, columns=cols)


def df_to_md(df):
    cols = list(df.columns)
    rows = ["| "+" | ".join(str(c) for c in cols)+" |",
            "| "+" | ".join("---" for _ in cols)+" |"]
    for _, row in df.iterrows():
        rows.append("| "+" | ".join("" if str(v)=="nan" else str(v) for v in row)+" |")
    return "\n".join(rows)


# ── Flowchart (Sankey) ────────────────────────────────────────────────────────

def make_flowchart(df_init, df_final, outpath):
    """4-stage Sankey: input → initial categorise → intermediate → final groups."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.path import Path
    except ImportError:
        print("matplotlib not available — flowchart skipped")
        return

    n = len(df_final)
    if n == 0:
        return

    def sc(df, prefix):
        return sum(v for k, v in df["STATUS"].value_counts().items() if k.startswith(prefix))

    # ── counts ──────────────────────────────────────────────────────────
    c_id    = sc(df_init, "ID_CONFIRMED")
    c_rid   = sc(df_init, "RESOLVED_BY_ID")
    c_ral   = sc(df_init, "RESOLVED_BY_ALIAS")
    c_uniq  = sc(df_init, "UNIQUE")
    c_ambig = sc(df_init, "AMBIGUOUS_UNRESOLVED")
    c_miss  = sc(df_init, "MISSING")

    c_clean = (sc(df_final, "ID_CONFIRMED") + sc(df_final, "RESOLVED_BY_ID") +
               sc(df_final, "RESOLVED_BY_ALIAS") + sc(df_final, "UNIQUE"))
    c_ckept = sc(df_final, "CONFLICT_KEPT")
    c_cdrop = sc(df_final, "CONFLICT_DROPPED")
    c_ambr  = sc(df_final, "AMBIGUOUS_UNRESOLVED")

    n_matched = c_clean + c_ckept
    n_dropped = c_cdrop + c_ambr
    n_nomatch = c_miss

    # ── colours ─────────────────────────────────────────────────────────
    C_GD = "#2e7d32"; C_GM = "#66bb6a"; C_GL = "#a5d6a7"
    C_YL = "#fff176"; C_SL = "#ffcdd2"; C_RL = "#e53935"
    C_RD = "#b71c1c"; C_GR = "#bdbdbd"
    # Dataset colours — avoid greens/reds already used for status
    DS_PALETTE = ["#4e79a7","#f28e2b","#76b7b2","#b07aa1",
                  "#edc948","#ff9da7","#9c755f","#bab0ac"]

    # ── Band definitions: (key, count, face, edge, label) BOTTOM → TOP ──
    # Ordering: NO MATCH at bottom, DROPPED above, MATCHED at top in all stages.
    # This keeps flow lines as uncrossed as possible.

    def _bands(*entries):
        return [(k, c, fc, ec, lb) for k, c, fc, ec, lb in entries if c > 0]

    # Stage 0: one band per dataset (sorted alphabetically, bottom → top)
    datasets   = sorted(df_init["DATASET"].unique())
    ds_counts  = {ds: (df_init["DATASET"] == ds).sum() for ds in datasets}
    ds_colors  = {ds: DS_PALETTE[i % len(DS_PALETTE)] for i, ds in enumerate(datasets)}
    s0 = _bands(*[
        (f"ds_{ds}", ds_counts[ds],
         ds_colors[ds], "#444444",
         f"{ds}\n({ds_counts[ds]:,})")
        for ds in datasets
    ])

    # Map s1 band keys → STATUS prefix for per-dataset flow computation
    S1_PREFIX = {"miss": "MISSING", "ambig": "AMBIGUOUS_UNRESOLVED",
                 "uniq": "UNIQUE",  "ral":   "RESOLVED_BY_ALIAS",
                 "rid":  "RESOLVED_BY_ID",   "id": "ID_CONFIRMED"}

    s1 = _bands(
        ("miss",  c_miss,  C_GR, "#757575", f"MISSING ({c_miss:,})"),
        ("ambig", c_ambig, C_YL, "#f9a825", f"AMBIG_UNRESOLVED\n({c_ambig:,})"),
        ("uniq",  c_uniq,  C_GL, C_GD,      f"UNIQUE ({c_uniq:,})"),
        ("ral",   c_ral,   C_GL, C_GD,      f"RESOLVED_BY_ALIAS\n({c_ral:,})"),
        ("rid",   c_rid,   C_GL, C_GD,      f"RESOLVED_BY_ID\n({c_rid:,})"),
        ("id",    c_id,    C_GL, C_GD,      f"ID_CONFIRMED\n({c_id:,})"),
    )

    s2 = _bands(
        ("miss2",  c_miss,  C_GR, "#757575", f"MISSING ({c_miss:,})"),
        ("ambr",   c_ambr,  C_RL, C_RD,      f"AMBIG_UNRESOLVED\n({c_ambr:,})"),
        ("cdrop",  c_cdrop, C_SL, C_RD,      f"CONFLICT_DROPPED\n({c_cdrop:,})"),
        ("ckept",  c_ckept, C_GM, C_GD,      f"CONFLICT_KEPT\n({c_ckept:,})"),
        ("clean",  c_clean, C_GL, C_GD,      f"direct resolved\n({c_clean:,})"),
    )

    s3 = _bands(
        ("nomatch", n_nomatch, C_RD, "#7f0000", f"NO MATCH\n{n_nomatch:,}\n{_pct(n_nomatch,n)}"),
        ("dropped", n_dropped, C_RL, C_RD,      f"DROPPED\n{n_dropped:,}\n{_pct(n_dropped,n)}"),
        ("matched", n_matched, C_GD, "#1b5e20",  f"MATCHED\n{n_matched:,}\n{_pct(n_matched,n)}"),
    )

    # ── Flow definitions: (src_key, dst_key, count, colour) ─────────────
    # Flows within each source band are ordered bottom-to-top matching the
    # destination band ordering — this minimises bezier crossings.

    # flows_01: each dataset band splits into s1 status bands.
    # Within each dataset band, flows are ordered bottom-to-top matching s1
    # so they fan out in the same direction and crossings are minimised.
    flows_01 = []
    for ds in datasets:
        df_ds  = df_init[df_init["DATASET"] == ds]
        ds_vc  = df_ds["STATUS"].value_counts()
        dc     = lambda pfx: sum(v for k,v in ds_vc.items() if k.startswith(pfx))
        for k, _cnt, fc, _ec, _lbl in s1:   # iterate s1 bottom→top
            cnt = dc(S1_PREFIX[k])
            if cnt:
                flows_01.append((f"ds_{ds}", k, cnt, fc))

    # Direct-resolved bands split proportionally into cdrop / ckept / clean.
    # Process s1 bands bottom-to-top (uniq → id), within each: cdrop first, clean last.
    n_dir = c_id + c_rid + c_ral + c_uniq
    def _dsplit(c):
        if n_dir == 0 or c == 0: return 0, 0, 0
        dr = round(c * c_cdrop / n_dir); kp = round(c * c_ckept / n_dir)
        return dr, kp, max(c - dr - kp, 0)

    flows_12 = []
    for k, cnt in [("uniq", c_uniq), ("ral", c_ral), ("rid", c_rid), ("id", c_id)]:
        dr, kp, cl = _dsplit(cnt)
        if dr: flows_12.append((k, "cdrop", dr, C_SL))
        if kp: flows_12.append((k, "ckept", kp, C_GM))
        if cl: flows_12.append((k, "clean", cl, C_GL))
    if c_ambr: flows_12.append(("ambig", "ambr",  c_ambr, C_RL))
    if c_miss: flows_12.append(("miss",  "miss2", c_miss, C_GR))

    flows_23 = []
    if c_miss:  flows_23.append(("miss2", "nomatch", c_miss,  C_RD))
    if c_ambr:  flows_23.append(("ambr",  "dropped", c_ambr,  C_RL))
    if c_cdrop: flows_23.append(("cdrop", "dropped", c_cdrop, C_RL))
    if c_ckept: flows_23.append(("ckept", "matched", c_ckept, C_GD))
    if c_clean: flows_23.append(("clean", "matched", c_clean, C_GD))

    # ── y-positions: stack bands bottom-to-top ───────────────────────────
    GAP = max(n * 0.015, 0.5)

    def _stack(bands):
        pos, y = {}, 0
        for i, (k, c, *_) in enumerate(bands):
            pos[k] = (y, y + c)
            y += c + (GAP if i < len(bands) - 1 else 0)
        return pos

    pos = [_stack(s) for s in (s0, s1, s2, s3)]
    ymax = max((p[max(p, key=lambda k: p[k][1])][1]) for p in pos if p)

    # ── Figure ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 7))
    bg = "#f7f7f7"
    fig.patch.set_facecolor(bg); ax.set_facecolor(bg)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.05 * ymax, 1.14 * ymax)
    ax.axis("off")

    XS = [0.06, 0.30, 0.67, 0.93]   # bar x-centres (0–1)
    BW = 0.022                        # bar half-width

    # ── Draw helpers ─────────────────────────────────────────────────────
    def draw_bars(xi, bands, band_pos, lbl_side):
        FONT, LS = 7.0, 1.3
        # Estimate one text line in data units so we can detect crowding.
        y0, y1  = ax.get_ylim()
        line_h  = (FONT * LS / 72) * ((y1 - y0) / fig.get_size_inches()[1])

        # Collect natural label positions (band centres) and per-label heights.
        items = []
        for k, cnt, fc, ec, lbl in bands:
            yb, yt = band_pos[k]
            items.append(dict(k=k, fc=fc, ec=ec, lbl=lbl, yb=yb, yt=yt,
                              ym=(yb + yt) / 2,
                              lbl_h=(lbl.count("\n") + 1) * line_h))

        # Bottom-to-top pass: push label positions apart when they would overlap,
        # requiring at least (sum of both half-heights + 20% buffer) between centres.
        placed = [it["ym"] for it in items]
        for i in range(1, len(items)):
            min_gap = (items[i-1]["lbl_h"] + items[i]["lbl_h"]) / 2 * 1.2
            if placed[i] < placed[i-1] + min_gap:
                placed[i] = placed[i-1] + min_gap

        # Draw bars then labels (with leader lines where staggered).
        for it, lbl_y in zip(items, placed):
            ax.add_patch(plt.Rectangle(
                (XS[xi] - BW, it["yb"]), 2*BW, it["yt"] - it["yb"],
                fc=it["fc"], ec=it["ec"], lw=1.2, zorder=3))

            if lbl_side == "center":
                ax.text(XS[xi], it["ym"], it["lbl"], ha="center", va="center",
                        fontsize=8.5, fontweight="bold", color="white",
                        linespacing=LS, zorder=4)
                continue

            sign  = 1 if lbl_side == "right" else -1
            x_edge = XS[xi] + sign * BW
            x_lbl  = x_edge + sign * 0.008
            ha     = "left" if lbl_side == "right" else "right"
            ax.text(x_lbl, lbl_y, it["lbl"], ha=ha, va="center",
                    fontsize=FONT, color="#1a1a1a", linespacing=LS, zorder=4)

            # Thin leader line from bar edge to staggered label position.
            if abs(lbl_y - it["ym"]) > line_h * 0.5:
                ax.plot([x_edge + sign * 0.002, x_edge + sign * 0.006],
                        [it["ym"], lbl_y],
                        color="#bbbbbb", lw=0.7, solid_capstyle="round", zorder=3)

    def bezier_band(x1, y1b, y1t, x2, y2b, y2t, color, alpha=0.40):
        cx = (x1 + x2) / 2
        verts = [(x1, y1t), (cx, y1t), (cx, y2t), (x2, y2t),
                 (x2, y2b), (cx, y2b), (cx, y1b), (x1, y1b), (x1, y1t)]
        codes = [Path.MOVETO,
                 Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO,
                 Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.CLOSEPOLY]
        ax.add_patch(mpatches.PathPatch(
            Path(verts, codes), fc=color, ec="none", alpha=alpha, zorder=2))

    def draw_flows(xi_src, src_pos, xi_dst, dst_pos, flows):
        src_off = {k: 0 for k in src_pos}
        dst_off = {k: 0 for k in dst_pos}
        for sk, dk, cnt, col in flows:
            if cnt <= 0 or sk not in src_pos or dk not in dst_pos:
                continue
            sy0 = src_pos[sk][0] + src_off[sk]
            dy0 = dst_pos[dk][0] + dst_off[dk]
            bezier_band(XS[xi_src] + BW, sy0, sy0 + cnt,
                        XS[xi_dst] - BW, dy0, dy0 + cnt, col)
            src_off[sk] += cnt
            dst_off[dk] += cnt

    # ── Render ────────────────────────────────────────────────────────────
    draw_bars(0, s0, pos[0], "left")
    draw_bars(1, s1, pos[1], "left")
    draw_bars(2, s2, pos[2], "right")
    draw_bars(3, s3, pos[3], "center")

    draw_flows(0, pos[0], 1, pos[1], flows_01)
    draw_flows(1, pos[1], 2, pos[2], flows_12)
    draw_flows(2, pos[2], 3, pos[3], flows_23)

    # ── Stage headers ─────────────────────────────────────────────────────
    for xi, title in enumerate(["INPUT\n(by dataset)", "INITIAL\nCATEGORISE",
                                 "INTERMEDIATE\n(conflict)", "FINAL"]):
        ax.text(XS[xi], ymax * 1.10, title, ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color="#333333")

    plt.tight_layout(pad=0.3)
    plt.savefig(outpath, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Flowchart saved → {outpath}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Resolve QRY→REF mapping from combined_summary.tsv.")
    p.add_argument("summary", nargs="?",
                   help="Input TSV (required unless --test)")
    p.add_argument("--aliases", default=DEFAULT_ALIASES,
                   help="Tab-delimited alias file (one group per line)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", default=None)
    p.add_argument("--test", action="store_true",
                   help="Run on built-in test dataset; no input file needed")
    args = p.parse_args()

    rng = random.Random(args.seed)

    if args.test:
        for f in (TEST_SUMMARY, TEST_ALIASES):
            if not f.exists():
                sys.exit(f"ERROR: test file not found: {f}")
        df = pd.read_csv(TEST_SUMMARY, sep="\t")
        ag = load_aliases(str(TEST_ALIASES))
        chart_path = str(TEST_DIR / "flowchart.png")
    else:
        if not args.summary:
            p.error("summary file required unless --test is given")
        stem       = os.path.splitext(args.summary)[0]
        outpath    = args.out or f"{stem}_resolved.tsv"
        stats_path = os.path.splitext(outpath)[0] + "_stats.tsv"
        chart_path = os.path.splitext(outpath)[0] + "_flowchart.png"
        df = pd.read_csv(args.summary, sep="\t")
        if miss := {"DATASET","QUERY","DUPLICATES"} - set(df.columns):
            sys.exit(f"ERROR: missing columns: {miss}")
        ag = load_aliases(args.aliases)

    # Collect known REF IDs before processing (for alias-safety validation)
    known_refs = _collect_ref_ids(df)

    # Run pipeline, preserving intermediate snapshots for the flowchart
    df_init = initial_categorise(df, ag)
    result  = check_surjectivity(df_init.copy(), rng)

    # ── Alias-safety validation ───────────────────────────────────────────
    sentinels = {"NA", "AMBIGUOUS"}
    bad = result[result["REF_MAPPED"].apply(
        lambda v: str(v) not in sentinels and str(v) not in known_refs)]
    if len(bad):
        print(f"\nWARNING: {len(bad)} rows have REF_MAPPED not in known REF IDs "
              f"(possible alias ID leakage):")
        print(bad[["QUERY","REF_MAPPED","STATUS"]].to_string(index=False))

    if args.test:
        W = 68

        # ── 1. INPUT ─────────────────────────────────────────────────────
        print("=" * W)
        print("TEST MODE — built-in dataset".center(W))
        print("=" * W)
        print(f"\n{'── INPUT ':─<{W}}")
        print(df.to_string(index=False))

        # ── 2. ALIASES ───────────────────────────────────────────────────
        print(f"\n{'── ALIASES ':─<{W}}")
        seen_groups = set()
        for g in ag.values():
            if g not in seen_groups:
                seen_groups.add(g)
                print("  " + "  ↔  ".join(sorted(g)))
        if not seen_groups:
            print("  (none loaded)")

        # ── 3. OUTPUT ────────────────────────────────────────────────────
        print(f"\n{'── OUTPUT ':─<{W}}")
        # Merge input DUPLICATES column alongside output columns for easy comparison
        out = (df[["DATASET","QUERY","DUPLICATES"]]
               .merge(result[["DATASET","QUERY","REF_MAPPED","STATUS","ALIAS_NOTE"]],
                      on=["DATASET","QUERY"], how="left"))
        # Suppress empty ALIAS_NOTE to keep the table readable
        out["ALIAS_NOTE"] = out["ALIAS_NOTE"].replace("", "—").fillna("—")
        print(out.to_string(index=False))

        # ── 4. CHECKS ────────────────────────────────────────────────────
        print(f"\n{'── CHECKS ':─<{W}}")
        expected = {
            "REF001":     ("ID_CONFIRMED",            "REF001"),
            "QRY001":     ("UNIQUE",                  "REF002"),
            "QRY_ALIAS":  ("RESOLVED_BY_ALIAS",       "REF_ALIAS_A"),
            "REF_T1":     ("RESOLVED_BY_ID",          "REF_T1"),
            "QRY_T_ALIAS":("RESOLVED_BY_ALIAS",       "REF_T2"),
            "QRY003":     ("AMBIGUOUS_UNRESOLVED",     "AMBIGUOUS"),
            "QRY004":     ("UNIQUE",                  "REF003"),
            "QRY005":     ("AMBIGUOUS_UNRESOLVED",    "AMBIGUOUS"),
            "QRY006":     ("UNIQUE",                  "REF005"),
            "QRY007":     ("UNIQUE",                  "REF006"),
            "QRY008":     ("MISSING",                 "NA"),
            "QRY011":     ("AMBIGUOUS_UNRESOLVED",    "AMBIGUOUS"),
        }
        out_map = result.set_index("QUERY")[["STATUS","REF_MAPPED"]]
        ok = fail = 0
        for qry, (exp_status, exp_ref) in expected.items():
            row = out_map.loc[qry] if qry in out_map.index else None
            got_status = row["STATUS"] if row is not None else "???"
            got_ref    = row["REF_MAPPED"] if row is not None else "???"
            mark = "OK " if got_status.startswith(exp_status) and got_ref == exp_ref else "FAIL"
            if mark == "OK ": ok += 1
            else:             fail += 1
            print(f"  {mark}  {qry:<14} status={got_status:<32} ref={got_ref}")

        # QRY009/QRY010: seed-dependent which wins; verify the pair as a unit
        conflict_pair = ["QRY009", "QRY010"]
        pair_statuses = [out_map.loc[q]["STATUS"] for q in conflict_pair if q in out_map.index]
        pair_refs     = [out_map.loc[q]["REF_MAPPED"] for q in conflict_pair if q in out_map.index]
        kept    = sum(1 for s in pair_statuses if s.startswith("CONFLICT_KEPT"))
        dropped = sum(1 for s in pair_statuses if s.startswith("CONFLICT_DROPPED"))
        ref_ok  = any(r == "REF007" for r in pair_refs) and any(r == "NA" for r in pair_refs)
        pair_ok = (kept == 1 and dropped == 1 and ref_ok)
        mark = "OK " if pair_ok else "FAIL"
        if pair_ok: ok += 1
        else:       fail += 1
        detail = "yes" if pair_ok else f"kept={kept} dropped={dropped} refs={pair_refs}"
        print(f"  {mark}  QRY009+QRY010   1×CONFLICT_KEPT + 1×CONFLICT_DROPPED: {detail}")

        verdict = "all passed" if fail == 0 else f"{fail} FAILED"
        print(f"\n  {ok}/{ok+fail} checks passed — {verdict}")

        # ── Write test/output.md ──────────────────────────────────────────
        alias_lines = []
        seen_groups = set()
        for g in ag.values():
            if g not in seen_groups:
                seen_groups.add(g)
                alias_lines.append("| " + " ↔ ".join(f"`{x}`" for x in sorted(g)) + " |")

        check_rows = []
        for qry, (exp_status, exp_ref) in expected.items():
            row = out_map.loc[qry] if qry in out_map.index else None
            got = row["STATUS"] if row is not None else "???"
            ref = row["REF_MAPPED"] if row is not None else "???"
            check_rows.append({"QUERY": qry, "STATUS": got, "REF_MAPPED": ref,
                                "": "✓" if got.startswith(exp_status) and ref == exp_ref else "✗"})
        check_rows.append({"QUERY": "QRY009+QRY010",
                           "STATUS": "1×CONFLICT_KEPT + 1×CONFLICT_DROPPED",
                           "REF_MAPPED": "—", "": "✓" if pair_ok else "✗"})
        checks_df = pd.DataFrame(check_rows)[["QUERY","STATUS","REF_MAPPED",""]]

        md = (f"### Input\n\n{df_to_md(df)}\n\n"
              f"### Aliases\n\n| Group |\n| --- |\n" +
              ("\n".join(alias_lines) if alias_lines else "| _(none)_ |") +
              f"\n\n### Output\n\n{df_to_md(out)}\n\n"
              f"### Checks\n\n{df_to_md(checks_df)}\n\n"
              f"**{ok}/{ok+fail} checks — {verdict}**\n")
        (TEST_DIR / "output.md").write_text(md)
        print(f"Test output written → {TEST_DIR / 'output.md'}")
    else:
        # ── Normal: write files and print summary ─────────────────────────
        result[["QUERY","REF_MAPPED","DATASET","STATUS","CANDIDATES","ALIAS_NOTE"]].to_csv(
            outpath, sep="\t", index=False)
        stats = build_stats_table(result)
        stats.to_csv(stats_path, sep="\t", index=False)

        n = len(result); counts = result["STATUS"].value_counts()
        def sc(pfx): return sum(v for k,v in counts.items() if k.startswith(pfx))
        print(f"\nInput: {args.summary}  →  {outpath}\n")
        for gh, _, statuses in SUMMARY_GROUPS:
            gn = sum(sc(s) for s,_ in statuses)
            print(f"  {gh:<10} {gn:>7,}  ({_pct(gn,n)})")
            for s,_ in statuses:
                cnt = sc(s)
                if cnt: print(f"    {s:<32} {cnt:>7,}  ({_pct(cnt,n)})")
        print(f"\n  {'TOTAL':<10} {n:>7,}")
        totals    = stats[stats["SECTION"]=="TOTALS"].drop(columns=["SECTION","STATUS"]).reset_index(drop=True)
        breakdown = stats[stats["SECTION"]=="BREAKDOWN"].drop(columns="SECTION").reset_index(drop=True)
        md_body   = f"## Mapping Totals\n\n{df_to_md(totals)}\n\n## Mapping Breakdown\n\n{df_to_md(breakdown)}\n"
        md_path   = os.path.splitext(stats_path)[0] + ".md"
        with open(md_path, "w") as fh:
            fh.write(md_body)
        print(md_body)

    make_flowchart(df_init, result, chart_path)


if __name__ == "__main__":
    main()
