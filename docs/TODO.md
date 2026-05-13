# TODO — Pancreatic Tumor ABM

Goal: Full replication of Akman Yıldız et al. (2021), then spatial and treatment extensions.

---

## Done

- [x] Architecture refactor (ParamGroup/JSON, AlwaysCopyToNew, thread-safe GlobalCensus, PopulationLogger)
- [x] ODE fidelity — exact Eqs. 2.1–2.6, Table 1 parameters, bilinear terms (no spurious Sat)
- [x] Global-mode validation — all 6 populations R² > 0.97 vs paper Fig. 2 (S=1e4, seed 42)
- [x] Run archive system — `scripts/save_run.py`, `runs/index.csv`, `scripts/run_experiment.sh`
- [x] Publication plots — `data-export/create-plots.py` (2×3 comparison grid, ABM + ODE + paper dots)
- [x] ODE reference solver — `scripts/ode_reference.py` (paper mode + abm mode)
- [x] Local-mode spatial extension:
  - [x] Gaussian random walk for immune cells (σ = 10 µm/step)
  - [x] Treg suppression via global R (avoids R_local=0 singularity)
  - [x] Carrying capacity via global C/P (avoids Poisson-pocket overshoot)
  - [x] Uniform initialization (sphere packing incompatible with density compensation)
  - [x] Local-mode v2 validated: all R² > 0.89 (days 7–82)

---

## Next — Drug Treatment (Paper Section 5)

Priority: **HIGH** — direct paper replication milestone.

The paper studies three drugs:
- **Abraxane** (nab-Paclitaxel): reduces PSC activation (lowers `p_boost_from_C` / `c_boost_from_P`)
- **Gemcitabine**: increases tumor death rate (raises `c_base_death` term or reduces `c_base_div`)
- **Anti-CD47**: blocks "don't eat me" signal → increases NK killing (`c_kill_by_N` ×10–100)

### Implementation steps

- [ ] Add treatment schedule to `SimParam`: `std::map<int, double> treatment_days` per drug type
- [ ] Add `ScheduledTreatmentOp` class (one per drug) — pattern: CARTopiaX `SpawnCart` in `src/utils/utils_aux.h` / `utils_aux.cc`. Operator checks `sim_day` against schedule and multiplies relevant rate params for that day.
- [ ] Register ops in `Simulate()` after scheduler setup
- [ ] Add treatment params to `params.json`: `abraxane_days`, `gem_days`, `anti_cd47_days`, `abraxane_dose`, etc.
- [ ] Run and validate against paper Fig. 3 (single drug) and Fig. 4 (combination)
- [ ] Archive treatment runs with descriptive notes

---

## Next — Cancer Stem Cell Extension (Paper Section 6)

Priority: **MEDIUM** — model extension beyond baseline.

Paper adds 7th variable S (CSCs):
```
dS/dt = k_s S (1 − a_s S) + φ C − d_s S
```
CSCs are drug-resistant and responsible for post-treatment relapse.

### Implementation steps

- [ ] New agent class `CancerStemCell` (`BDM_AGENT_HEADER`, same pattern as `TumorCell`)
- [ ] New behavior `CSCBehavior`:
  - Division: `k_s * (1 - S/K_S)` — higher rate than PCC, smaller carrying cap
  - Death: `d_s` — much lower than PCC
  - Drug resistance: drug treatment ops skip or reduce kill probability on CSC
- [ ] Conversion terms:
  - PCC → CSC: `TumorBehavior` converts at rate φ per step (call `this->ConvertTo<CancerStemCell>()` or remove + spawn)
  - CSC → PCC: `CSCBehavior` differentiates at rate `d_phi` per step
- [ ] Add CSC to `GlobalCensus` (S column), `SourceBehavior`, `PopulationLogger`
- [ ] Add CSC params to `SimParam` and `params.json`
- [ ] Validate against paper Fig. 5 (relapse after Gemcitabine treatment)

---

## Next — PRCC Sensitivity Analysis (Paper Section 4)

Priority: **MEDIUM** — paper methodology replication.

Paper reports PRCC of C and E at days 7, 35, 70, 150 across 1000 LHS samples.

### Implementation steps

- [ ] `scripts/lhs_sample.py` — Latin Hypercube Sampling over all Table 1 params in [0.5×, 2×]:
  - Output: `lhs/params_000.json` through `lhs/params_999.json`
  - Use `scipy.stats.qmc.LatinHypercube`
- [ ] `scripts/run_lhs.sh` — parallel batch runner:
  ```bash
  for i in $(seq -w 0 999); do
    ./build/pancreatic_tumor_new --config lhs/params_$i.json \
      --output lhs/out_$i.csv &
  done
  wait
  ```
  Requires adding `--config` and `--output` CLI flags to `main.cc`
- [ ] `scripts/prcc.py` — PRCC computation:
  - Load all `lhs/out_*.csv` + `lhs/params_*.json` into a DataFrame
  - Rank-transform inputs and outputs
  - Partial correlation controlling for all other params
  - Plot PRCC bar charts vs paper Fig. 3
- [ ] Validate: top sensitive params for C should be `a_c`/`k_c`/`μ_c`; for E should be `r₁`/`δ_e`

---

## Code Quality (Low Priority)

- [ ] **PopulationLogger flush** — add `csv_.flush()` after each `WriteRow()` (KI-07: prevents data loss on crash)
- [ ] **`std::endl` in `ReportPopCounts`** — replace `"\n"` so piped runs show progress in real time (stdout buffering)
- [ ] **Replace dynamic_cast census** with `TimeSeries` + `Counter<T>` for performance at >50k agents (KI-05)
- [ ] **Sphere initialization for local mode** — `CreateSphereOfTumorCells()` HCP pattern (CARTopiaX `src/utils/utils_aux.cc`); requires spatial calibration of density compensation

---

## Reference: Paper Sections → Implementation Map

| Paper section | Content | Status |
|--------------|---------|--------|
| Sec. 2 — ODE model | Eqs. 2.1–2.6, Table 1 | ✅ Complete |
| Sec. 3 — Baseline dynamics | Fig. 2, untreated 100 days | ✅ Complete |
| Sec. 4 — Sensitivity | PRCC, Fig. 3 | ☐ Pending |
| Sec. 5 — Drug treatment | Abraxane, Gemcitabine, Anti-CD47, Figs. 3–4 | ☐ Pending |
| Sec. 6 — CSC extension | 7th variable S, relapse dynamics, Fig. 5 | ☐ Pending |
