#!/usr/bin/env python3
"""
post_treatment_suite.py — aggregate a treatment production suite and export
paper-Fig-5-style comparison plots + a replication-quality metrics table.

Consumes the layout produced by run_treatment_production.py:
    runs/treatment_prod/<suite>/<protocol>/<TS>_<protocol>_s<seed>/
        populations.csv, ode_reference.csv, fit_metrics_summary.csv, params.json

For each protocol it computes ABM mean±std across seeds and compares against the
treatment ODE (paper Section 5, Eqs 5.1-5.7 — the digitized paper Fig. 5 curves
are NOT available, so the ODE is the reference target).

Outputs (in <suite>/analysis/ by default):
    treatment_fig5_tumor.png      one panel per protocol: tumor C(t), ABM vs ODE,
                                  treatment windows shaded  (paper Fig. 5 style)
    treatment_compare_tumor.png   all protocols' tumor C(t) on one axis (efficacy)
    treatment_<proto>_allpops.png per-protocol 2x3 grid of all 6 populations
    treatment_drugs.png           drug concentrations M_gem/M_abr(t) per protocol
    treatment_metrics_summary.csv per protocol x population: R2/MAPE mean±std,
                                  endpoint counts  (ABM vs ODE)

Usage:
    python3 scripts/post_treatment_suite.py --suite-dir runs/treatment_prod/<suite>
    python3 scripts/post_treatment_suite.py --suite-dir <suite> --protocols gem,abr
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Population display metadata (matches create-plots.py / post_validation_suite.py)
POPS = ["c", "p", "e", "n", "h", "r"]
POP_LABEL = {"c": "Tumor (C)", "p": "PSC (P)", "e": "CD8 (E)",
             "n": "NK (N)", "h": "Helper (H)", "r": "Treg (R)"}
POP_COLOR = {"c": "#c0392b", "p": "#e67e22", "e": "#2980b9",
             "n": "#27ae60", "h": "#8e44ad", "r": "#7f8c8d"}
# Distinct colors for protocol overlays in the cross-protocol comparison
PROTO_COLOR = {
    "untreated": "#2c3e50", "acd47": "#16a085", "gem": "#2980b9",
    "abr": "#e67e22", "abr_acd47": "#8e44ad", "abr_csc": "#c0392b",
}
# Drug-window shading colors
DRUG_SHADE = {"gem": ("#2980b9", "Gemcitabine"),
              "abr": ("#e67e22", "Abraxane"),
              "acd47": ("#16a085", "Anti-CD47")}


def discover_protocols(suite_dir: Path):
    """Return {protocol: [run_dir, ...]} for every protocol subfolder with runs."""
    protos = {}
    for pdir in sorted(suite_dir.iterdir()):
        if not pdir.is_dir() or pdir.name == "analysis":
            continue
        runs = sorted([d for d in pdir.iterdir()
                       if d.is_dir() and (d / "populations.csv").exists()])
        if runs:
            protos[pdir.name] = runs
    return protos


def aggregate(run_dirs):
    """Load all seeds for one protocol; return mean/std/min/max + ODE + params.

    Only fully-complete runs are aggregated: a still-running seed produces a
    short populations.csv, and including it would truncate every seed to its
    length. Keep runs at the maximal (full) length and drop shorter ones.
    """
    loaded = []  # (run_dir, df)
    for rd in run_dirs:
        df = pd.read_csv(rd / "populations.csv")
        df.columns = df.columns.str.lower().str.strip()
        loaded.append((rd, df))
    full_len = max(len(df) for _, df in loaded)
    complete = [(rd, df) for rd, df in loaded if len(df) == full_len]
    dropped = len(loaded) - len(complete)
    if dropped:
        print(f"    (skipped {dropped} incomplete run(s); using {len(complete)} seeds)")
    run_dirs = [rd for rd, _ in complete]
    dfs = [df for _, df in complete]

    days = dfs[0]["days"].values
    stack = {p: np.vstack([df[p].values for df in dfs]) for p in POPS}
    agg = {
        "days": days,
        "n_seeds": len(dfs),
        "mean": {p: stack[p].mean(axis=0) for p in POPS},
        "std":  {p: stack[p].std(axis=0) for p in POPS},
        "min":  {p: stack[p].min(axis=0) for p in POPS},
        "max":  {p: stack[p].max(axis=0) for p in POPS},
    }
    # Drug concentrations (mean is enough — deterministic schedule)
    for m in ("m_gem", "m_abr"):
        if m in dfs[0].columns:
            agg[m] = np.vstack([df[m].values for df in dfs]).mean(axis=0)

    # ODE reference + params from the first run
    ode_csv = run_dirs[0] / "ode_reference.csv"
    if ode_csv.exists():
        ode = pd.read_csv(ode_csv)
        ode.columns = ode.columns.str.lower().str.strip()
        agg["ode"] = ode
    params = json.load(open(run_dirs[0] / "params.json"))
    agg["params"] = params
    agg["csc"] = bool(params.get("csc_enable"))
    return agg


def shade_windows(ax, params, label=False):
    """Shade [treat_start_day, <drug>_end_day] for each active drug."""
    start = params.get("treat_start_day", 14.0)
    active = []
    if params.get("treat_gem"):   active.append(("gem", params.get("gem_end_day", 56.0)))
    if params.get("treat_abr"):   active.append(("abr", params.get("abr_end_day", 28.0)))
    if params.get("treat_acd47"): active.append(("acd47", params.get("acd47_end_day", 35.0)))
    for drug, end in active:
        color, name = DRUG_SHADE[drug]
        ax.axvspan(start, end, color=color, alpha=0.10,
                   label=(f"{name} window" if label else None))


def _finalize_log_axis(ax, values):
    vals = [v for v in values if v is not None and np.isfinite(v) and v > 0]
    if vals:
        ax.set_ylim(max(0.5, min(vals) * 0.5), max(vals) * 2.0)
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)


# --- Fig. 5 style: tumor C(t) per protocol -------------------------------
def plot_fig5_tumor(aggs, out_path):
    protos = list(aggs.keys())
    n = len(protos)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 4.2 * nrow), squeeze=False)
    for i, proto in enumerate(protos):
        ax = axes[i // ncol][i % ncol]
        a = aggs[proto]
        days = a["days"]
        allv = []
        # ABM mean +/- std
        ax.fill_between(days, np.maximum(a["mean"]["c"] - a["std"]["c"], 1e-9),
                        a["mean"]["c"] + a["std"]["c"],
                        color=POP_COLOR["c"], alpha=0.25, zorder=2,
                        label=f"ABM mean±std (n={a['n_seeds']})")
        ax.plot(days, a["mean"]["c"], color=POP_COLOR["c"], lw=1.6, ls="--",
                zorder=3, label="ABM mean")
        allv += list(a["mean"]["c"])
        # ODE reference
        if "ode" in a:
            ax.plot(a["ode"]["days"], a["ode"]["c"], color="black", lw=1.4,
                    ls=":", zorder=4, label="ODE reference (Eqs 5.1-5.7)")
            allv += list(a["ode"]["c"])
        shade_windows(ax, a["params"], label=True)
        title = proto + ("  [CSC relapse]" if a["csc"] else "")
        if a["csc"]:
            title += "\n(ODE ref lacks CSC — qualitative only)"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Day (paper)"); ax.set_ylabel("Tumor cells C")
        _finalize_log_axis(ax, allv)
        ax.legend(fontsize=7, loc="best")
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Treatment response — tumor C(t): ABM vs ODE reference (Fig. 5 style)",
                 fontsize=13, y=1.0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path.name}")


# --- Cross-protocol tumor efficacy comparison ----------------------------
def plot_compare_tumor(aggs, out_path):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    allv = []
    for proto, a in aggs.items():
        color = PROTO_COLOR.get(proto, None)
        ax.plot(a["days"], a["mean"]["c"], color=color, lw=1.8,
                label=f"{proto} (n={a['n_seeds']})")
        allv += list(a["mean"]["c"])
    # Shade EACH drug's dosing window (union across all protocols), once per drug,
    # with its own colour + a labelled legend entry and an end-of-dosing marker.
    windows = {}  # drug -> (start, end)
    for a in aggs.values():
        p = a["params"]
        start = p.get("treat_start_day", 14.0)
        if p.get("treat_gem"):   windows.setdefault("gem",   (start, p.get("gem_end_day", 56.0)))
        if p.get("treat_abr"):   windows.setdefault("abr",   (start, p.get("abr_end_day", 28.0)))
        if p.get("treat_acd47"): windows.setdefault("acd47", (start, p.get("acd47_end_day", 35.0)))
    # Draw widest first so the (nested) narrower windows stay visible on top.
    for drug, (start, end) in sorted(windows.items(), key=lambda kv: -(kv[1][1] - kv[1][0])):
        color, name = DRUG_SHADE[drug]
        ax.axvspan(start, end, color=color, alpha=0.13, zorder=0,
                   label=f"{name} dosing (d{start:.0f}–{end:.0f})")
        ax.axvline(end, color=color, lw=1.1, ls="--", alpha=0.6, zorder=0)

    ax.set_title("Tumor response by protocol — ABM mean C(t)\n"
                 "shaded = each drug's dosing window (dashed = dosing ends)")
    ax.set_xlabel("Day (paper)"); ax.set_ylabel("Tumor cells C")
    _finalize_log_axis(ax, allv)
    ax.legend(fontsize=8, loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path.name}")


# --- Cross-protocol tumor comparison WITH the ODE / math model overlaid ---
def plot_compare_tumor_ode(aggs, out_path):
    """Same as plot_compare_tumor (all protocols' tumor C(t) on one axis, with
    each drug's dosing window shaded), but ALSO overlays each protocol's ODE /
    mathematical-model solution as a dotted line in the same colour."""
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(12, 7))
    allv = []
    for proto, a in aggs.items():
        color = PROTO_COLOR.get(proto, None)
        # ABM mean — solid continuous line
        ax.plot(a["days"], a["mean"]["c"], color=color, lw=1.8,
                label=f"{proto} (n={a['n_seeds']})", zorder=3)
        allv += list(a["mean"]["c"])
        # ODE / mathematical model — dotted line, same colour
        if "ode" in a and "c" in a["ode"].columns:
            ax.plot(a["ode"]["days"], a["ode"]["c"], color=color, lw=2.0,
                    ls=":", alpha=0.9, zorder=4)
            allv += list(a["ode"]["c"])

    # Shade each drug's dosing window (union across protocols), once per drug.
    windows = {}
    for a in aggs.values():
        p = a["params"]; start = p.get("treat_start_day", 14.0)
        if p.get("treat_gem"):   windows.setdefault("gem",   (start, p.get("gem_end_day", 56.0)))
        if p.get("treat_abr"):   windows.setdefault("abr",   (start, p.get("abr_end_day", 28.0)))
        if p.get("treat_acd47"): windows.setdefault("acd47", (start, p.get("acd47_end_day", 35.0)))
    for drug, (start, end) in sorted(windows.items(), key=lambda kv: -(kv[1][1] - kv[1][0])):
        wc, name = DRUG_SHADE[drug]
        ax.axvspan(start, end, color=wc, alpha=0.13, zorder=0,
                   label=f"{name} dosing (d{start:.0f}–{end:.0f})")
        ax.axvline(end, color=wc, lw=1.1, ls="--", alpha=0.6, zorder=0)

    ax.set_title("Tumor response by protocol — ABM (solid) vs mathematical model / ODE (dotted)\n"
                 "shaded = each drug's dosing window")
    ax.set_xlabel("Day (paper)"); ax.set_ylabel("Tumor cells C")
    _finalize_log_axis(ax, allv)
    # Legend 1: protocol colours + drug windows
    leg1 = ax.legend(fontsize=8, loc="upper left", ncol=2)
    ax.add_artist(leg1)
    # Legend 2: line style -> which curve
    style = [Line2D([0], [0], color="0.25", lw=1.8, ls="-", label="ABM (mean of seeds)"),
             Line2D([0], [0], color="0.25", lw=2.0, ls=":", label="ODE / math model (theory)")]
    ax.legend(handles=style, fontsize=8, loc="lower right", title="Line style")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path.name}")


# --- Per-protocol all-population grid ------------------------------------
def plot_allpops(proto, a, out_path):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for idx, p in enumerate(POPS):
        ax = axes[idx // 3][idx % 3]
        days = a["days"]
        allv = []
        ax.fill_between(days, np.maximum(a["mean"][p] - a["std"][p], 1e-9),
                        a["mean"][p] + a["std"][p],
                        color=POP_COLOR[p], alpha=0.25, zorder=2, label="ABM mean±std")
        ax.plot(days, a["mean"][p], color=POP_COLOR[p], lw=1.4, ls="--",
                zorder=3, label="ABM mean")
        allv += list(a["mean"][p])
        if "ode" in a and p in a["ode"].columns:
            ax.plot(a["ode"]["days"], a["ode"][p], color="black", lw=1.2, ls=":",
                    zorder=4, label="ODE")
            allv += list(a["ode"][p])
        shade_windows(ax, a["params"])
        ax.set_title(POP_LABEL[p]); ax.set_xlabel("Day"); ax.set_ylabel("cells")
        _finalize_log_axis(ax, allv)
        ax.legend(fontsize=7)
    fig.suptitle(f"Protocol '{proto}' — ABM (n={a['n_seeds']}) vs ODE, all populations",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path.name}")


# --- Drug concentration confirmation -------------------------------------
def plot_drugs(aggs, out_path):
    have = {p: a for p, a in aggs.items() if ("m_gem" in a or "m_abr" in a)
            and any(a["params"].get(k) for k in ("treat_gem", "treat_abr"))}
    if not have:
        return
    fig, ax = plt.subplots(figsize=(11, 6))
    for proto, a in have.items():
        if a["params"].get("treat_gem") and "m_gem" in a:
            ax.plot(a["days"], a["m_gem"], lw=1.4, label=f"{proto}: M_gem")
        if a["params"].get("treat_abr") and "m_abr" in a:
            ax.plot(a["days"], a["m_abr"], lw=1.4, ls="--", label=f"{proto}: M_abr")
    ax.set_title("Drug concentrations M(t) — PK schedule confirmation")
    ax.set_xlabel("Day (paper)"); ax.set_ylabel("Concentration M")
    ax.grid(True, alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path.name}")


# --- Replication metrics table (aggregate per-seed fit_metrics) -----------
def build_metrics(suite_dir, protos):
    rows = []
    for proto, run_dirs in protos.items():
        per_pop = {}   # population -> {"R2":[...], "MAPE":[...]}
        endpoint_c = []
        csc = False
        n_used = 0
        for rd in run_dirs:
            fm = rd / "fit_metrics_summary.csv"
            params = rd / "params.json"
            # Skip still-running seeds: save_run.py writes these only at completion.
            if not fm.exists() or not params.exists():
                continue
            n_used += 1
            m = pd.read_csv(fm)
            for _, r in m.iterrows():
                d = per_pop.setdefault(r["population"], {"R2": [], "MAPE": []})
                d["R2"].append(r["R2"]); d["MAPE"].append(r["MAPE_%"])
            pop = pd.read_csv(rd / "populations.csv")
            pop.columns = pop.columns.str.lower().str.strip()
            endpoint_c.append(pop["c"].values[-1])
            csc = csc or bool(json.load(open(params)).get("csc_enable"))
        for population, d in per_pop.items():
            rows.append({
                "protocol": proto,
                "population": population,
                "n_seeds": n_used,
                "R2_mean": round(float(np.mean(d["R2"])), 4),
                "R2_std": round(float(np.std(d["R2"])), 4),
                "MAPE_mean": round(float(np.mean(d["MAPE"])), 2),
                "endpoint_C_mean": round(float(np.mean(endpoint_c)), 1),
                "endpoint_C_std": round(float(np.std(endpoint_c)), 1),
                "ode_ref_valid": (not csc),  # abr_csc ODE lacks CSC relapse term
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite-dir", required=True)
    ap.add_argument("--out-dir", default=None,
                    help="Default: <suite-dir>/analysis")
    ap.add_argument("--protocols", default=None,
                    help="Comma list to restrict (default: all found)")
    args = ap.parse_args()

    suite_dir = Path(args.suite_dir).resolve()
    out_dir = Path(args.out_dir) if args.out_dir else suite_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    protos = discover_protocols(suite_dir)
    if args.protocols:
        keep = {p.strip() for p in args.protocols.split(",")}
        protos = {k: v for k, v in protos.items() if k in keep}
    if not protos:
        print(f"ERROR: no protocol runs found under {suite_dir}")
        return 1

    print(f"Suite: {suite_dir}")
    print(f"Protocols: {[(p, len(r)) for p, r in protos.items()]}")

    aggs = {p: aggregate(r) for p, r in protos.items()}

    print("\nGenerating plots...")
    plot_fig5_tumor(aggs, out_dir / "treatment_fig5_tumor.png")
    plot_compare_tumor(aggs, out_dir / "treatment_compare_tumor.png")
    plot_compare_tumor_ode(aggs, out_dir / "treatment_compare_tumor_ode.png")
    plot_drugs(aggs, out_dir / "treatment_drugs.png")
    for proto, a in aggs.items():
        plot_allpops(proto, a, out_dir / f"treatment_{proto}_allpops.png")

    print("\nComputing replication metrics (ABM vs ODE)...")
    metrics = build_metrics(suite_dir, protos)
    metrics_csv = out_dir / "treatment_metrics_summary.csv"
    metrics.to_csv(metrics_csv, index=False)

    # Print a compact tumor-focused replication table
    tumor = metrics[metrics["population"] == "Tumor"].copy()
    print("\n  Tumor replication (ABM vs paper ODE):")
    print("  " + "-" * 68)
    print(f"  {'protocol':12s} {'n':>3s} {'R2_mean':>9s} {'R2_std':>8s} "
          f"{'end_C':>10s} {'ode_ref':>8s}")
    for _, r in tumor.iterrows():
        flag = "ok" if r["ode_ref_valid"] else "NO-CSC"
        print(f"  {r['protocol']:12s} {int(r['n_seeds']):>3d} "
              f"{r['R2_mean']:>9.3f} {r['R2_std']:>8.3f} "
              f"{r['endpoint_C_mean']:>10.0f} {flag:>8s}")
    print("  " + "-" * 68)
    print(f"\nWrote metrics -> {metrics_csv}")
    print(f"Plots + metrics in: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
