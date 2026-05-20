# SynthDiD.jl

`SynthDiD.jl` is a Julia implementation of synthetic difference-in-differences
from Arkhangelsky, Athey, Hirshberg, Imbens, and Wager (2021). It follows the
basic interface of the R [`synthdid`](https://github.com/synth-inference/synthdid)
package.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/xao/SynthDiD.jl")  # or Pkg.develop(path="...")
```

## Quick Start

```julia
using SynthDiD

# Load California Proposition 99 data
df = california_prop99()
setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)

# Three estimators
τ_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0)  # -15.60
τ_sc   = sc_estimate(setup.Y, setup.N0, setup.T0)        # -19.62
τ_did  = did_estimate(setup.Y, setup.N0, setup.T0)       # -27.35

# Standard errors
se_val = se(τ_sdid; method=:placebo, replications=200)

# Effect curve (per-period treatment effects)
curve = effect_curve(τ_sdid)
```

## API

### Estimators

- `synthdid_estimate(Y, N0, T0; X=nothing, units=nothing, times=nothing, ...)` — Synthetic Diff-in-Diff (Algorithm 1)
- `sc_estimate(Y, N0, T0)` — Synthetic Control
- `did_estimate(Y, N0, T0)` — Standard Diff-in-Diff

All three return a `SynthDiDEstimate` with fields: `estimate`, `estimator`, `weights`, `setup`, `opts`.

### Data Preparation

- `panel_matrices(df, unit, time, outcome, treatment; covariates=Symbol[])` — convert long panel to matrix format
- `california_prop99()` — load bundled dataset

### Inference

- `vcov(est; method, replications)` — variance estimate (`:bootstrap`, `:jackknife`, `:placebo`)
- `se(est; method, replications)` — standard error

### Diagnostics

- `effect_curve(est)` — per-period treatment effects
- `placebo(est)` — placebo estimate using pre-treatment data
- `synthdid_controls(est; mass=0.9, weight_type="omega")` — top control weights as DataFrame

### Plotting

With `using Plots`:
- `plot(est)` — parallel trends with time weight bars (supports `overlay` parameter)
- `plot([est1, est2, est3])` — side-by-side comparison

### Covariate Adjustment

Pass covariates through `panel_matrices` and the `X` keyword:

```julia
setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated;
                       covariates=[:income, :population])
τ = synthdid_estimate(setup.Y, setup.N0, setup.T0; X=setup.X)
```

## Method

The estimator solves a weighted TWFE regression:

$$\hat\tau^{sdid} = \arg\min_{\tau,\mu,\alpha,\beta} \sum_{i,t} (Y_{it} - \mu - \alpha_i - \beta_t - W_{it}\tau)^2 \hat\omega_i \hat\lambda_t$$

where $\hat\omega$ (unit weights) and $\hat\lambda$ (time weights) are estimated via Frank-Wolfe optimization on the unit simplex with ridge regularization.

The useful way to remember the estimator is that SC, DiD, and SDiD differ only
in the unit and time weights:
- **SC**: $\omega$ optimized, $\lambda = 0$, no unit intercept
- **DiD**: $\omega = 1/N_0$ (uniform), $\lambda = 1/T_0$ (uniform)
- **SDiD**: Both $\omega$ and $\lambda$ optimized

## Vignettes

- [Introduction](https://xiangao.github.io/SynthDiD.jl/tutorials/01_introduction/) — California Prop 99 walkthrough
- [More on Plotting](https://xiangao.github.io/SynthDiD.jl/tutorials/02_more_plotting/) — Plot customization and weight tables
- [Paper Results](https://xiangao.github.io/SynthDiD.jl/tutorials/03_paper_results/) — Replicating paper results

Full documentation: **https://xiangao.github.io/SynthDiD.jl/**

Local render: `julia --project=docs docs/make.jl`

## Reference

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). Synthetic Difference-in-Differences. *American Economic Review*, 111(12), 4088-4118.

## Dependencies

- DataFrames.jl, CSV.jl, Distributions.jl, RecipesBase.jl
- Statistics, LinearAlgebra, Random (stdlib)
