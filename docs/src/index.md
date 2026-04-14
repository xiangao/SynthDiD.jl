# SynthDiD.jl

`SynthDiD.jl` is a Julia implementation of **Synthetic Difference-in-Differences** (Arkhangelsky, Athey, Hirshberg, Imbens, and Wager, 2021). It provides a compact interface for synthetic DiD, synthetic control, and standard DiD estimation from balanced panel data.

## Features

- **Three estimators**: `synthdid_estimate`, `sc_estimate`, and `did_estimate`
- **Panel reshaping utilities**: `panel_matrices` converts long data into the matrix representation used by the estimators
- **Inference tools**: `vcov`, `se`, and placebo-based diagnostics
- **Plotting support**: built-in plot recipes plus accessors for custom visualizations

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/xao/SynthDiD.jl")
```

## Quick Start

```julia
using SynthDiD

df = california_prop99()
setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)

tau_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0)
tau_sc = sc_estimate(setup.Y, setup.N0, setup.T0)
tau_did = did_estimate(setup.Y, setup.N0, setup.T0)
```

## Tutorials

- [Introduction](tutorials/01_introduction.md): California Proposition 99 walkthrough and estimator overview
- [More on Plotting](tutorials/02_more_plotting.md): plot customization, weight summaries, and accessors
- [Paper Results](tutorials/03_paper_results.md): replication-style estimates and placebo simulations
