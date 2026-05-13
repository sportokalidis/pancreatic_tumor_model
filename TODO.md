# Implementation Plan

Goal: ODE-replicating ABM at conference quality. All 6 populations should track
the reference paper (Akman Yıldız 2021) with R² > 0.7 and MAPE < 50%.

---

## In Progress

*(none)*

---

## To Do

### Medium priority (improve accuracy / future extension)

- [ ] **Replace dynamic_cast census with TimeSeries + Counter<T>** — remove
  the `ForEachAgent` scan, use BioDynaMo's built-in population tracking

- [ ] **Sphere initialization for tumor** — start C cells in a compact sphere
  at the domain center (groundwork for spatial extension)

- [ ] **Drug treatment scenarios** (paper Sec. 5) — Abraxane, Gemcitabine,
  Anti-CD47; add treatment parameters and time-dependent kill terms

---

## Done

- [x] **Disable mechanical forces** — `UnscheduleOp("mechanical forces")` in
  `Simulate()`; eliminates domain overcrowding artifact (KI-02)

- [x] **Add constant immune source terms** (KI-01) — separate per-cell
  proliferation from absolute influx (`a_e`, `a_n`, `a_h`, `a_r`) in behavior
  code; implement influx via `SourceBehavior` on reporter agent

- [x] **Fix P0** (KI-04) — `P0: 2 → 27` in params.json

- [x] **Rescale inactivation/suppression K constants** (KI-03) — 
  `e_inact_by_C_K`, `n_inact_by_C_K` → 100; `e_suppr_by_R_K`, `n_suppr_by_R_K` → 10;
  values now meaningful at ABM scale (hundreds of cells)

- [x] **Architecture refactor** — `ParamGroup`/JSON config, `AlwaysCopyToNew()`,
  `Initialize(NewAgentEvent)`, thread-safe `GlobalCensus`, `PopulationLogger`

- [x] **Analysis documentation** — `docs/analysis.md`, `docs/known_issues.md`,
  `docs/changelog.md`

- [x] **Automated validation** — `scripts/validate.py`, `scripts/run_and_validate.sh`

- [x] **Expanded test suite** — unit tests for `Sat`, `ProbFromRate`,
  `GlobalCensus`, `SimParam`, 2-day smoke test

- [x] **ODE reference Python solver** — `scripts/ode_reference.py` integrates
  Eqs. 2.1-2.6 with same params; outputs `data-export/ode_reference.csv` and
  `data-export/ode_vs_abm.png`; ABM tracks ODE within 1–15% at all timepoints

- [x] **Publication-quality comparison plots** — `data-export/create-plots.py`
  generates `comparison_grid.png` (2×3 panel) and `comparison_combined.png`
  (overview) overlaying ABM, ODE, and paper reference on the same axes

- [x] **Paper ODE corrected** — `scripts/ode_reference.py` now implements exact
  Eqs. 2.1-2.6 with Table 1 parameters scaled to ABM (S=1e5). `--mode paper`
  regenerates `*_scaled_global.csv` so reference dots and ODE line agree.
  Dynamics: E peaks ~17k at day 37 then collapses; NK→4; H→14; R→818;
  C and P reach carrying capacity. Key scaling: r1→34500, g_e/g_n/g_h/g_r→674490.

- [x] **Run archive & index system** — `scripts/save_run.py` archives each run
  into `runs/YYYYMMDD_HHMMSS_s{seed}/` (params, output, ODE ref, metrics, plots,
  run_info.json); `runs/index.csv` accumulates one row per run with all params +
  metrics for sensitivity analysis; `scripts/run_experiment.sh` is the full
  build→run→archive pipeline

- [x] **Align ABM to Table 1 — conference quality achieved** — params.json and
  C++ behaviors updated to exact Eqs. 2.1-2.6 (bilinear terms replacing Hill
  approximations; r1=34500, lambda_p≈0, g=674490, all correct scaling).
  Result: C R²=0.9999 MAPE=0.6%, P R²=0.995 MAPE=18%, E R²=0.999 MAPE=19%,
  N R²=0.997 MAPE=13%, H R²=0.988 MAPE=15%, R R²=0.975 MAPE=7%.
  All populations exceed R²>0.97, well above the R²>0.7 conference target.
