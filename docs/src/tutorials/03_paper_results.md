# Paper Results

```@meta
CurrentModule = SynthDiD
```

This tutorial reproduces the package's main Proposition 99 estimates and a simple placebo exercise in the style of the original Synthetic Difference-in-Differences paper.

## Setup

```@example paper
using SynthDiD
using Plots
using DataFrames
using Statistics
using Random
using Printf

Random.seed!(12345)
```

## California Proposition 99 Application

Estimate the three models on the bundled panel:

```@example paper
df = california_prop99()
setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
years = collect(setup.times)

tau_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0;
                             units=setup.units, times=years)
tau_sc = sc_estimate(setup.Y, setup.N0, setup.T0;
                     units=setup.units, times=years)
tau_did = did_estimate(setup.Y, setup.N0, setup.T0;
                       units=setup.units, times=years)

estimates = [tau_did, tau_sc, tau_sdid]
est_names = ["DiD", "SC", "SDiD"]
```

Summarize point estimates and placebo standard errors:

```@example paper
point_est = [e.estimate for e in estimates]

Random.seed!(12345)
se_placebo = [se(e; method=:placebo, replications=200) for e in estimates]

println("\n" * "=" ^ 60)
@printf("%-20s %12s %12s\n", "Estimator", "tau", "Placebo SE")
println("=" ^ 60)
for (name, est, se_val) in zip(est_names, point_est, se_placebo)
    @printf("%-20s %12.2f %12.2f\n", name, est, se_val)
end
println("=" ^ 60)
```

Plot the three fitted paths together:

```@example paper
p1 = plot(estimates; size=(1000, 350))
savefig(p1, "paper_estimates.svg")
nothing
```

![](paper_estimates.svg)

Inspect the largest control and time weights:

```@example paper
unit_weights = synthdid_controls(estimates; weight_type="omega", mass=0.95)
time_weights = synthdid_controls(estimates; weight_type="lambda", mass=0.95)

println(unit_weights)
println(time_weights)
```

## Placebo Simulations

A basic placebo routine can be used to summarize how noisy the estimator is under resampling from the donor pool:

```@example paper
function run_placebo_sims(Y, N0, T0; n_sims=100)
    N1 = size(Y, 1) - N0
    estimates = Float64[]

    for sim in 1:n_sims
        ind = rand(1:N0, N0)
        n0_new = length(ind) - N1
        if n0_new < 1
            continue
        end

        try
            est = synthdid_estimate(Y[ind, :], n0_new, T0)
            push!(estimates, est.estimate)
        catch
            continue
        end
    end

    estimates
end

Random.seed!(12345)
placebo_est = run_placebo_sims(setup.Y, setup.N0, setup.T0; n_sims=100)
```

Summarize and visualize the placebo distribution:

```@example paper
@printf("Mean: %1.2f\n", mean(placebo_est))
@printf("Std:  %1.2f\n", std(placebo_est))
@printf("RMSE: %1.2f\n", sqrt(mean(placebo_est .^ 2)))
```

```@example paper
p2 = histogram(placebo_est;
    bins=20,
    label="Placebo estimates",
    xlabel="Estimated treatment effect",
    ylabel="Frequency",
    title="Placebo Simulation (100 replications)",
    color=:steelblue,
    alpha=0.7)

vline!([tau_sdid.estimate];
    label="SDiD estimate",
    color=:firebrick,
    linewidth=2)
savefig(p2, "paper_placebo_histogram.svg")
nothing
```

![](paper_placebo_histogram.svg)

## Interpretation

Across these examples, SDiD produces a more conservative estimate than synthetic control and a much smaller effect than plain DiD, reflecting the balance obtained by jointly choosing unit and time weights.
