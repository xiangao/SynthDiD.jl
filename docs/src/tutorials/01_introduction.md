# Introduction to SynthDiD.jl

`SynthDiD.jl` implements three related estimators for policy evaluation with panel data:

- **Synthetic Diff-in-Diff (SDiD)** optimizes both unit weights `ω` and time weights `λ`
- **Synthetic Control (SC)** optimizes unit weights only
- **Diff-in-Diff (DiD)** uses uniform weights across controls and pre-treatment periods

This tutorial uses the California Proposition 99 tobacco-tax example bundled with the package.

## Setup

```julia
using SynthDiD
using Plots
using DataFrames
using Random
using Printf

Random.seed!(12345)
```

## California Proposition 99

Load the long-format panel and inspect the first few rows:

```julia
df = california_prop99()
first(df, 5)
```

Convert the panel to the matrix representation expected by the estimators:

```julia
setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
years = collect(setup.times)

println("Y: $(size(setup.Y)), N0 = $(setup.N0), T0 = $(setup.T0)")
```

## Point Estimates

Estimate all three models on the same treated-control setup:

```julia
tau_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0;
                             units=setup.units, times=years)
tau_sc = sc_estimate(setup.Y, setup.N0, setup.T0;
                     units=setup.units, times=years)
tau_did = did_estimate(setup.Y, setup.N0, setup.T0;
                       units=setup.units, times=years)

println(tau_sdid)
println(tau_sc)
println(tau_did)
```

Compare the resulting treatment effects directly:

```julia
println("-" ^ 50)
@printf("%-30s %10s\n", "Estimator", "tau")
println("-" ^ 50)
for (name, est) in [("Diff-in-Diff", tau_did),
                    ("Synthetic Control", tau_sc),
                    ("Synthetic Diff-in-Diff", tau_sdid)]
    @printf("%-30s %10.2f\n", name, est.estimate)
end
println("-" ^ 50)
```

## Standard Errors

With one treated unit, placebo inference is the appropriate variance estimator:

```julia
Random.seed!(12345)
se_val = se(tau_sdid; method=:placebo, replications=200)

@printf("Point estimate: %1.2f\n", tau_sdid.estimate)
@printf("Placebo SE:     %1.2f\n", se_val)
@printf("95%% CI:         (%1.2f, %1.2f)\n",
        tau_sdid.estimate - 1.96 * se_val,
        tau_sdid.estimate + 1.96 * se_val)
```

## Plotting

The built-in recipe plots treated and synthetic trajectories together:

```julia
plot(tau_sdid)
```

You can also compare all three estimators side by side:

```julia
plot([tau_did, tau_sc, tau_sdid])
```

## Effect Curves and Weights

Inspect the post-treatment effect path:

```julia
curve = effect_curve(tau_sdid)
post_years = years[(setup.T0 + 1):end]

bar(post_years, curve; label="Per-period effect")
```

Summarize the most influential control units:

```julia
omega = tau_sdid.weights.omega
control_names = setup.units[1:setup.N0]
top_idx = sortperm(omega, rev=true)[1:10]

for i in top_idx
    @printf("%-20s %0.4f\n", control_names[i], omega[i])
end
```
