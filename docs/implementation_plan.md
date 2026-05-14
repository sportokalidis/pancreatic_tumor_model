# Implementation Plan — Pancreatic Tumor ABM

Phased roadmap from current state to a complete, validated BioDynaMo model.

**Status (2026-05-14):**
- Phase 1 (Architecture): ✅ Complete
- Phase 2 (ODE Fidelity): ✅ Complete — global mode all R² > 0.97; local mode v2 all R² > 0.89
- Phase 3 (Spatial Init): 🔶 Partial — uniform init working; sphere init deferred (dc incompatibility)
- Phase 4 (Validation): ✅ Complete — automated validate.py, run archive, comparison plots
- Phase 5 (Treatment): 🔶 Implemented (branch: drug-treatment) — needs dose calibration + ODE reference for validation
- Phase 6 (CSC): ☐ Pending
- Phase 7 (Sensitivity): ☐ Pending

---

## Phase 1 — Code Architecture Refactoring ✅ COMPLETE

The current code has all agents, behaviors, parameters and simulation logic in a single `.h` file using fragile patterns. This must be fixed before extending the model.

### 1.1 Adopt BioDynaMo `ParamGroup` for parameters

**Current problem:** `struct Params` with a global singleton `inline Params* P()`.
**Target pattern:** CARTopiaX `SimParam : public ParamGroup` with JSON loading.

Steps:
- Create `src/params/sim_param.h` and `src/params/sim_param.cc`
- Inherit from `bdm::ParamGroup`, use `BDM_PARAM_GROUP_HEADER`
- Move all parameters from `Params` struct into `SimParam`
- Add `LoadParams(filename)` that reads `params.json` (no recompile needed)
- Add `PrintParams()` for debugging
- Create `params.json` with all default values
- Reference in code via `Simulation::GetActive()->GetParam()->Get<SimParam>()`

**Why critical:** Without JSON config, every parameter change requires a full recompile. The singleton `P()` is not thread-safe and breaks parallel runs.

### 1.2 Split into modular files

```
src/
  agents/
    tumor_cell.h / tumor_cell.cc          — TumorCell class
    stellate_cell.h / stellate_cell.cc    — StellateCell
    effector_tcell.h / effector_tcell.cc  — EffectorTCell
    nk_cell.h / nk_cell.cc               — NKCell
    helper_tcell.h / helper_tcell.cc      — HelperTCell
    treg_cell.h / treg_cell.cc            — TRegCell
  behaviors/
    tumor_behavior.h / .cc
    psc_behavior.h / .cc
    effector_behavior.h / .cc
    nk_behavior.h / .cc
    helper_behavior.h / .cc
    treg_behavior.h / .cc
  utils/
    census.h / census.cc                  — GlobalCensus + LocalCounter
    reporter.h / reporter.cc              — CSV output
    init_populations.h / .cc             — Population seeding
  params/
    sim_param.h / sim_param.cc
  pancreatic_tumor_model.h               — Simulate() only
  pancreatic_tumor_model.cc             — main()
```

### 1.3 Fix cell division pattern

**Current problem:** After `c->Divide()`, code does `dynamic_cast<TumorCell*>(d)` and manually adds behavior. This is fragile — if BioDynaMo creates the daughter differently, the cast returns null.

**Target pattern (CARTopiaX):**
```cpp
// In behavior constructor:
TumorBehavior() { AlwaysCopyToNew(); }  // automatically copied to daughter

// In TumorCell::Initialize(const NewAgentEvent& event):
// Copy parent fields to daughter cell
void Initialize(const NewAgentEvent& event) override {
  Base::Initialize(event);
  // inherit color, state, etc. from mother
}
```

Steps:
- Add `AlwaysCopyToNew()` call in every behavior constructor
- Implement `Initialize(const NewAgentEvent& event)` in each agent class
- Remove manual `dynamic_cast` + `AddBehavior` after `Divide()` calls

### 1.4 Fix static file stream in reporter

**Current problem:** `static std::ofstream csv("data-export/populations.csv")` inside `ReportPopCounts::Run()` — hardcoded path, only one instance possible, not closeable.

**Target:** Move CSV management to a dedicated `Reporter` class initialized in `Simulate()` with a configurable output path (from params).

---

## Phase 2 — ODE Fidelity ✅ COMPLETE

Align ABM rules exactly with the paper ODEs before adding spatial features.

### 2.1 Audit and fix each behavior against the ODE

For each cell type, verify:
1. Division/birth rate formula matches ODE growth term
2. Death rate formula matches ODE death term
3. All interaction terms from the ODE are present
4. No extra terms not in the ODE (e.g., the `gate_C_K` gating in global mode is not in the paper)

**Key discrepancies to fix in global mode:**

| ODE term | Current code | Issue |
|----------|-------------|-------|
| `d_c EC/(1+r_1 R)` | `c_kill_by_E * Sat(E, K) / (1 + c_R_blocks_E * R)` | `Sat(E,K)` instead of direct E. In ODE, it's linear in E. |
| `b_c NC` | `c_kill_by_N * Sat(N, K)` | Same: linear in ODE, saturating in ABM |
| `r_e NC` (E recruitment by NK-tumor contact) | Missing | Not implemented |
| cytokine term `f̂_e (β_E E + β_N N + β_H H) E` | Missing | Not implemented |
| `r_e NC` for NK recruitment | Missing | Not implemented |
| Treg basal source `a` | `r_base_src` | Check sign and magnitude match |
| E basal: `a_e` (constant recruitment) | `e_base_birth * crowdE` | Should be constant source, not crowding-gated |

### 2.2 Parameter calibration

Map the paper's Table 1 parameters to ABM rates:
- Compute scaled rates: `rate_ABM = rate_paper * (ABM_pop / paper_pop)`
- Verify that the ABM steady state (after ~100 days) matches paper Fig. 2

### 2.3 Validate against paper Fig. 2

Target: reproduce the qualitative shape of all 6 population curves from Fig. 2:
- C and P rise to ~10⁹ carrying capacity
- E, N, H decline over time (immune suppression by tumor)
- R gradually rises

Run 5 independent seeds and compare mean trajectories.

---

## Phase 3 — Spatial Initialization 🔶 PARTIAL

The paper assumes a homogeneous spatial distribution, but BioDynaMo enables 3D structure.

### 3.1 Tumor sphere initialization

Instead of placing all C cells randomly in the full 3D box:
- Initialize C cells in a spherical cluster of radius `r_init` (e.g., 60–100 µm)
- Place immune cells (E, N, H, R) at the periphery or uniformly in the rest of the box
- Place PSC cells (P) near the tumor cluster

Reference: `CreateSphereOfTumorCells()` in CARTopiaX `utils/utils_aux.cc`.

### 3.2 Simulation domain scaling

Current domain: 150 × 150 × 150 µm with 6 µm cell radius.
Maximum packing: ~(150/12)³ ≈ 1953 cells — roughly correct for ~2000 initial agents.

Review and document the domain size choice for each parameter scenario.

### 3.3 Neighbor search radius

Ensure BioDynaMo's `interaction_length` (neighbor search radius) covers `local_radius_um` correctly:
```cpp
param->interaction_length = 2 * P()->local_radius_um + 1;
```
This is currently commented out in `Simulate()`.

---

## Phase 4 — Validation & Comparison ✅ COMPLETE

### 4.1 Error metric

Implement normalized RMSE between ABM population trajectories and paper ODE reference:
```
RMSE(X) = sqrt(mean((X_ABM(t) - X_ODE(t))² / X_ODE(t)²))
```
The `data-export/calc-error.py` script exists — verify it uses this metric.

### 4.2 Multiple runs

Run N=10 independent seeds and compute:
- Mean trajectory per population
- 95% confidence interval band
- RMSE vs paper reference for each population

### 4.3 Population scaling

The ABM uses ~2000 agents vs paper's ~10⁸ cells. Document and justify the scaling factor used. Ensure qualitative behaviors (not absolute values) match.

---

## Phase 5 — Treatment Modeling ☐ NEXT MILESTONE (Priority: HIGH)

Replicate paper Sections 4–5 (chemotherapy and immunotherapy).

**Reference pattern:** CARTopiaX `SpawnCart` (`src/utils/utils_aux.h` / `utils_aux.cc`) — `std::map<int, int>` day→dose schedule in SimParam; `operator()()` checks current simulation day and applies effect.

### 5.1 Drug injection mechanism

Implement a `ScheduledTreatmentOp` operation:
- Configurable treatment schedule: `map<int day, double dose>` per drug type in SimParam
- Three drug effects (modify rate params for duration of treatment):
  - **Abraxane** (nab-Paclitaxel): reduces PSC activation (`p_boost_from_C`, `c_boost_from_P`)
  - **Gemcitabine**: increases tumor death rate (adds kill term proportional to dividing cells)
  - **Anti-CD47**: increases NK killing of tumor (`c_kill_by_N` multiplied by dose factor)

### 5.2 Drug PK/PD model

Add drug concentration dynamics:
- Drug injected on specific days (bolus or continuous)
- Exponential decay: `D(t) = D₀ * exp(-λ_D * t)`
- Effect on relevant parameters as a function of D(t)

### 5.3 Treatment scenarios

Reproduce paper's finding: drug dose > drug frequency for efficacy.

---

## Phase 6 — Cancer Stem Cell Extension (Priority: MEDIUM-LOW)

Implement paper Section 6: add CSC (S) as a 7th variable.

### CSC ODE (from paper)

```
dS/dt = k_s S (1 − a_s S) + φ C − d_s S
```

Where:
- k_s: CSC self-renewal rate
- a_s: inverse carrying capacity
- φ C: differentiation from PCC (bidirectional)
- d_s: CSC death

### ABM implementation

- New agent class `CancerStemCell`
- Behavior: higher proliferation rate, drug resistance (reduced kill probability)
- PCC → CSC conversion at low rate φ
- CSC → PCC differentiation

---

## Phase 7 — Sensitivity Analysis (Priority: MEDIUM)

Replicate paper Section 4 (PRCC sensitivity analysis).

### 7.1 LHS sampling script

Create a Python script that:
- Defines parameter ranges [0.5×, 2×] for each parameter
- Generates N=1000 Latin Hypercube Samples
- Writes one `params.json` per sample

### 7.2 Parallel run infrastructure

```bash
#!/bin/bash
for i in $(seq 0 999); do
  ./build/pancreatic_tumor_model --config lhs/params_$i.json &
done
wait
```

### 7.3 PRCC computation

Python script to:
- Read all `output_i.csv` files
- Compute PRCC for C and E at t = 7, 35, 70, 150 days
- Reproduce paper Fig. 3

---

## Implementation Priority Order

```
Phase 1 (Architecture)  ──► ✅ Done
Phase 2 (ODE Fidelity)  ──► ✅ Done
Phase 4 (Validation)    ──► ✅ Done
Phase 3 (Spatial)       ──► 🔶 Partial (uniform init done; sphere deferred)
Phase 5 (Treatment)     ──► ☐ NEXT — paper Sec. 5
Phase 6 (CSC)           ──► ☐ paper Sec. 6
Phase 7 (Sensitivity)   ──► ☐ paper Sec. 4
```

---

## Milestones

| Milestone | Criteria | Status |
|-----------|----------|--------|
| M1: Clean architecture | JSON params, modular files, correct Divide pattern | ✅ |
| M2: ODE replication | All 6 populations R² > 0.97 vs paper Fig. 2 | ✅ (global: all > 0.97; local v2: all > 0.89) |
| M3: Spatial model | Tumor sphere init, local-mode validated | 🔶 Local validated; sphere init deferred |
| M4: Treatment | Abraxane+Gemcitabine+Anti-CD47 reproduce paper Sec. 5 | 🔶 Implemented; validation pending |
| M5: CSC | Relapse dynamics with CSCs match paper Sec. 6 | ☐ |
| M6: Sensitivity | PRCC outputs match paper Fig. 3 | ☐ |
