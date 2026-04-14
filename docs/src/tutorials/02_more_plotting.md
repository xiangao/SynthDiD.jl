# More on Plotting

`SynthDiD.jl` stores unit labels, time labels, and weights alongside each estimate, which makes it straightforward to go beyond the default plot recipe.

## Setup

```julia
using SynthDiD
using Plots
using DataFrames
using Random

Random.seed!(12345)
```

## Real Labels for Units and Time

Pass `units` and `times` when estimating so the plots use the original panel labels:

```julia
df = california_prop99()
setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
years = collect(setup.times)

tau_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0;
                             units=setup.units, times=years)
tau_sc = sc_estimate(setup.Y, setup.N0, setup.T0;
                     units=setup.units, times=years)
tau_did = did_estimate(setup.Y, setup.N0, setup.T0;
                       units=setup.units, times=years)
```

## Built-in Plot Recipe

The quickest way to visualize the treatment effect is the package recipe:

```julia
plot(tau_sdid)
```

To compare estimators, pass a vector:

```julia
plot([tau_did, tau_sc, tau_sdid])
```

## Custom Plot Construction

Because the estimate object exposes `weights` and `setup`, custom plots are easy to build:

```julia
Y = setup.Y
N0 = setup.N0
T0 = setup.T0
omega = tau_sdid.weights.omega

synth = vec(omega' * Y[1:N0, :])
california = vec(Y[N0 + 1, :])

p = plot(years, synth;
    label="Synthetic Control",
    color=:steelblue,
    linewidth=2,
    xlabel="Year",
    ylabel="Packs per Capita")

plot!(p, years, california;
    label="California",
    color=:firebrick,
    linewidth=2)

vline!(p, [years[T0] + 0.5];
    label="",
    color=:gray40,
    linestyle=:dash)
```

## Weight Summaries

The `synthdid_controls` helper returns the most important unit weights as a `DataFrame`:

```julia
ctrl = synthdid_controls([tau_sdid, tau_sc, tau_did]; mass=0.9)
println(ctrl)
```

You can also inspect time weights directly:

```julia
lambda = tau_sdid.weights.lambda
bar(years[1:T0], lambda; label="Time weights")
```
