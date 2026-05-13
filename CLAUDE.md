# Pancreatic Tumor Model — CLAUDE.md

Agent-based replication of the mathematical pancreatic cancer model by Akman Yıldız et al. (2021), built on the BioDynaMo platform.

---

## Project Goal

Replicate the ODE dynamics (temporal evolution of 6 cell populations) of the reference paper using a 3D agent-based model in BioDynaMo. Then extend with:
- Spatial structure (tumor sphere initialization)
- Drug treatment scenarios (Abraxane, Gemcitabine, Anti-CD47)
- Cancer Stem Cell (CSC) hypothesis (paper Section 6)
- Sensitivity analysis (PRCC)

---

## Reference Papers

| Role | Reference |
|------|-----------|
| Primary model | Akman Yıldız T, Köse E, Elliott SL. *Mathematical modeling of pancreatic cancer treatment with cancer stem cells.* J Biol Syst. 2021;29(4):799–832. DOI: 10.1142/S0218339021500182. PDF: `pancreatic-tumor-mathematical-model-paper.pdf` |
| BioDynaMo platform | Breitwieser L et al. *BioDynaMo: a modular platform for high-performance agent-based simulation.* Bioinformatics. 2022;38(2):453–460. |
| Closest reference ABM | CARTopiaX (`/home/stavport/Documents/dev/bdm-projects/CARTopiaX`) — CAR T-cell therapy ABM on BioDynaMo, same pattern to follow |

---

## Cell Populations (Model Variables)

| Symbol | Cell Type | BioDynaMo Agent |
|--------|-----------|-----------------|
| C | Pancreatic Cancer Cells (PCC) | `TumorCell` |
| P | Pancreatic Stellate Cells (PSC) | `StellateCell` |
| E | CD8⁺ T cells (Effectors) | `EffectorTCell` |
| N | Natural Killer cells (NK) | `NKCell` |
| H | Helper T cells (CD4⁺) | `HelperTCell` |
| R | Regulatory T cells (Tregs) | `TRegCell` |

---

## Key Source Files

| File | Purpose |
|------|---------|
| `src/pancreatic_tumor_model.h` | All agents, behaviors, parameters, census, reporter, `Simulate()` |
| `src/pancreatic_tumor_model.cc` | `main()` entry point |
| `data-export/populations.csv` | Simulation output (step, days, C, P, E, N, H, R, total) |
| `data-export/create-plots.py` | Python plot script |
| `data-export/calc-error.py` | Error metric computation against reference data |
| `docs/model_equations.md` | ODE equations from the reference paper |
| `docs/implementation_plan.md` | Phased development roadmap |
| `docs/suggestions.md` | Code robustness improvements needed |

---

## Reference ABM Projects (BioDynaMo)

| Path | What to learn |
|------|--------------|
| `/home/stavport/Documents/dev/bdm-projects/CARTopiaX/` | Closest analogue: modular file structure, `ParamGroup`+JSON config, `AlwaysCopyToNew()`, diffusion grids, custom forces, `Initialize(NewAgentEvent)` pattern |
| `/home/stavport/Documents/dev/bdm-projects/lung-project/` | Population reporter via `custom_ops`, `sim-param.h` pattern |
| `/home/stavport/Documents/dev/bdm-projects/bdm-paper-examples/` | BioDynaMo API usage examples |

---

## BioDynaMo Source (for API reference)

- Core classes: `/home/stavport/Documents/dev/biodynamo/src/core/`
- Demos: `/home/stavport/Documents/dev/biodynamo/demo/`

---

## Build & Run

```bash
# Source BioDynaMo environment first
source /path/to/biodynamo/build/bin/thisbdm.sh

# Build
biodynamo build
# or manually:
mkdir build && cd build && cmake .. && make -j$(nproc)

# Run
biodynamo run
# or:
./build/pancreatic_tumor_model

# Output written to: data-export/populations.csv
```

---

## Model Parameters (from paper Table 1, scaled for ABM)

See `docs/model_equations.md` for the full parameter table.

### Initial Conditions (ABM-scaled, current code)

| Variable | Paper value | ABM value | Scale factor |
|----------|------------|-----------|-------------|
| C₀ | 4.886 × 10⁷ | 474 | ~10⁵ |
| P₀ | 2.7362 × 10⁵ | 2 | ~10³ |
| E₀ | 4.2684 × 10⁷ | 431 | ~10⁵ |
| N₀ | 2.3531 × 10⁷ | 189 | ~10⁵ |
| H₀ | 1.0343 × 10⁸ | 839 | ~10⁵ |
| R₀ | 7.7570 × 10⁶ | 65 | ~10⁵ |

### Time mapping

- Simulation time unit: **days** (paper)
- BioDynaMo step unit: **minutes** (`dt_minutes = 1.0`)
- Steps per day: 1440
- `ProbFromRate(rate_per_day, dt_day)` converts ODE rates to per-step probabilities via `1 - exp(-rate * dt)`

---

## Interaction Mode

Two modes controlled by `Params::use_local_counts`:

- **Global mode** (`false`): counts all agents across simulation → matches ODE mean-field assumption. Preferred for ODE replication.
- **Local mode** (`true`): counts neighbors in `local_radius_um` sphere → adds spatial heterogeneity. Requires `RecomputeCapsForLocal()`.

Use global mode first to validate against paper Fig. 2, then explore local mode.

---

## Known Issues & TODOs

See `docs/suggestions.md` for a prioritized list of robustness improvements needed.

---

## Current Branch Status

- Branch `local-Neighbohood-counter`: experimenting with local neighborhood counting
- Branch `master`: global-mode baseline
- Output data in `data-export/` validated against paper reference curves (qualitative match)
