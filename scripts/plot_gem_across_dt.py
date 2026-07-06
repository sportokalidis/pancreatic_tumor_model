#!/usr/bin/env python3
"""Side-by-side gemcitabine C(t) across dt — illustrates the PK-resolution effect.
Shows ABM (each dt) vs the dt-independent ODE, with the gem treatment window shaded.
"""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs" / "treatment_prod"
SUITES = {
    "dt=24h": "20260705_204550_S1e4_E10x_dt24h_dosefix",
    "dt=6h":  "20260705_205646_S1e4_E10x_dt6h_dosefix",
    "dt=1h":  "20260705_205646_S1e4_E10x_dt1h_dosefix",
}


def load(path):
    r = list(csv.DictReader(open(path)))
    return ([float(x["days"]) for x in r], [float(x["C"]) for x in r])


fig, ax = plt.subplots(figsize=(11, 6.5))
colors = {"dt=24h": "#c0392b", "dt=6h": "#e67e22", "dt=1h": "#27ae60"}
ode_plotted = False
for label, suite in SUITES.items():
    gem = list((RUNS / suite / "gem").glob("*/populations.csv"))
    if not gem:
        continue
    d, c = load(gem[0])
    ax.plot(d, c, color=colors[label], lw=1.8, label=f"ABM gem {label}")
    if not ode_plotted:
        ode = list((RUNS / suite / "gem").glob("*/ode_reference.csv"))
        if ode:
            od, oc = load(ode[0])
            ax.plot(od, oc, color="black", lw=1.6, ls=":", label="ODE (paper Eqs 5.1-5.7)")
            ode_plotted = True

ax.axvspan(14, 56, color="#2980b9", alpha=0.08, label="Gemcitabine window")
ax.set_yscale("log")
ax.set_xlabel("Day (paper)")
ax.set_ylabel("Tumor cells C")
ax.set_title("Gemcitabine (dose=0.1) — tumor response vs timestep\n"
             "Coarse dt over-kills the 3h-half-life drug; dt=1h resolves it and matches the ODE")
ax.grid(True, which="both", alpha=0.25)
ax.legend(fontsize=9, loc="lower right")
out = RUNS / "gem_across_dt.png"
fig.tight_layout()
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"wrote {out}")
