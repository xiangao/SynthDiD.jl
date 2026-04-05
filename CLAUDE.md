# SynthDiD.jl — AI Assistant Instructions

## Overview

Julia implementation of Synthetic Difference-in-Differences (Arkhangelsky et al. 2021), translated from the R `synthdid` package. Implements Algorithm 1 (SDiD), plus SC and DiD as special cases, with optional covariate adjustment.

## Project Structure

```
SynthDiD.jl/
├── Project.toml              # Package metadata and dependencies
├── src/SynthDiD.jl           # Single-file module (all code + plot recipes)
├── test/runtests.jl          # 91 tests: estimates, invariances, variance, diagnostics
├── data/california_prop99.csv # Bundled dataset (semicolon-delimited)
├── vignettes/
│   ├── Project.toml          # Separate env (Pkg.develop's parent package)
│   ├── introduction.qmd     # Main walkthrough (California Prop 99)
│   ├── more-plotting.qmd    # Plot customization and synthdid_controls
│   └── paper-results.qmd    # Paper replication (Tables, placebo sims)
├── output/                   # Rendered vignettes (HTML)
├── CLAUDE.md                 # This file
└── README.md
```

## Architecture

Single module file with these components:
- **Structs**: `SynthDiDEstimate`, `SynthDiDWeights`, `SynthDiDSetup`, `SynthDiDOpts`
  - Setup stores `Y`, `N0`, `T0`, `X` (covariates), `units`, `times` (labels)
  - Weights stores `omega`, `lambda`, `beta` (covariate coefficients)
- **Panel utilities**: `panel_matrices()` (with optional `covariates` kwarg), `collapsed_form()`
- **Solver**: `fw_step()` (Frank-Wolfe on simplex), `sc_weight_fw()`, `sc_weight_fw_covariates()` (alternating optimization for ω, λ, β)
- **Estimators**: `synthdid_estimate()`, `sc_estimate()`, `did_estimate()`
- **Inference**: `vcov()` with bootstrap/jackknife/placebo methods
- **Diagnostics**: `effect_curve()`, `placebo()`, `synthdid_controls()`, `california_prop99()`
- **Plot recipes** (RecipesBase): `@recipe` for single estimate (parallel trends + time weights, with optional `overlay`) and vector of estimates (side-by-side comparison). Uses Plots.jl, not Makie.

## R Source Reference

Translated from `~/projects/claude/repo_cloned/synthdid/R/`. Key mapping:
- `synthdid.R` → estimators
- `solver.R` → Frank-Wolfe optimizer (including `sc.weight.fw.covariates`)
- `utils.R` → `panel_matrices`, `collapsed_form`
- `vcov.R` → variance estimation (Algorithms 2-4)

## Testing

```bash
cd ~/projects/software/SynthDiD.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Reference values (California Prop 99): SDiD ≈ -15.60, SC ≈ -19.62, DiD ≈ -27.35

## Vignettes

Render with Quarto (uses separate vignette environment):
```bash
cd ~/projects/software/SynthDiD.jl/vignettes
JULIA_PROJECT=$(pwd) quarto render introduction.qmd
JULIA_PROJECT=$(pwd) quarto render more-plotting.qmd
JULIA_PROJECT=$(pwd) quarto render paper-results.qmd
```

## Not Implemented

- `synthdid_units_plot()` — R's per-unit contribution scatter plot (not ported due to requiring Plots.jl as module dependency)
