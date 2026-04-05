using Test
using SynthDiD
using Statistics
using Random
using Distributions
using LinearAlgebra
using DataFrames

@testset "SynthDiD.jl" begin

    @testset "California Prop 99" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)

        @test setup.N0 == 38
        @test setup.T0 == 19
        @test size(setup.Y) == (39, 31)

        τ_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0)
        τ_sc = sc_estimate(setup.Y, setup.N0, setup.T0)
        τ_did = did_estimate(setup.Y, setup.N0, setup.T0)

        # Known R values (approximate)
        @test isapprox(τ_sdid.estimate, -15.60, atol=1.0)
        @test isapprox(τ_sc.estimate, -19.62, atol=1.0)
        @test isapprox(τ_did.estimate, -27.35, atol=1.0)

        println("SDiD: ", round(τ_sdid.estimate, digits=2))
        println("SC:   ", round(τ_sc.estimate, digits=2))
        println("DiD:  ", round(τ_did.estimate, digits=2))
    end

    @testset "Effect curve" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ = synthdid_estimate(setup.Y, setup.N0, setup.T0)
        curve = effect_curve(τ)

        @test length(curve) == size(setup.Y, 2) - setup.T0
        # Average of effect curve should equal the point estimate
        T1 = length(curve)
        @test isapprox(mean(curve), τ.estimate, atol=0.01)
    end

    @testset "Placebo" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ = synthdid_estimate(setup.Y, setup.N0, setup.T0)
        τ_plac = placebo(τ)

        # Placebo estimate should be small (no treatment in pre-period)
        @test abs(τ_plac.estimate) < abs(τ.estimate)
    end

    # ─── Invariance tests (from R test suite) ─────────────────────────

    function random_low_rank(; rng=Random.default_rng())
        n0, n1, T0, T1 = 100, 10, 120, 20
        n = n0 + n1
        T = T0 + T1
        tau = 1.0
        sigma = 0.5
        rank = 2
        rho = 0.7

        var_mat = [rho^abs(i-j) for i in 1:T, j in 1:T]
        W = Float64.((1:n .> n0) * (1:T .> T0)')
        U = Float64.([rand(rng, Poisson(sqrt(s) / sqrt(n))) for s in sample(rng, 1:n, n), _ in 1:rank])
        V = Float64.([rand(rng, Poisson(sqrt(t) / sqrt(T))) for t in 1:T, _ in 1:rank])
        alpha = 10.0 .* sample(rng, 1:n, n) ./ n .* ones(1, T)
        beta = ones(n, 1) .* (10.0 .* (1:T)' ./ T)
        mu = U * V' + alpha + beta
        error_mat = rand(rng, MvNormal(zeros(T), var_mat), n)' .* sigma
        Y = mu + tau .* W + error_mat

        (Y=Y, L=mu, N0=n0, T0=T0)
    end

    @testset "Column fixed effect invariance" begin
        rng = MersenneTwister(42)
        setup = random_low_rank(rng=rng)
        bt = 2.0 .* ones(size(setup.Y, 1), 1) * (1:size(setup.Y, 2))'

        for estimator in [synthdid_estimate, sc_estimate, did_estimate]
            est1 = estimator(setup.Y, setup.N0, setup.T0)
            est2 = estimator(setup.Y + bt, setup.N0, setup.T0)
            @test isapprox(est1.estimate, est2.estimate, atol=1e-8)
        end
    end

    @testset "Row fixed effect invariance (not SC)" begin
        rng = MersenneTwister(42)
        setup = random_low_rank(rng=rng)
        ai = 2.5 .* (1:size(setup.Y, 1)) .* ones(1, size(setup.Y, 2))

        for estimator in [synthdid_estimate, did_estimate]
            est1 = estimator(setup.Y, setup.N0, setup.T0)
            est2 = estimator(setup.Y + ai, setup.N0, setup.T0)
            @test isapprox(est1.estimate, est2.estimate, atol=1e-8)
        end
    end

    @testset "Scaling invariance" begin
        rng = MersenneTwister(42)
        setup = random_low_rank(rng=rng)

        for c in [0.5, 2.0, 100.0]
            for estimator in [synthdid_estimate, sc_estimate, did_estimate]
                est1 = estimator(setup.Y, setup.N0, setup.T0)
                est2 = estimator(c .* setup.Y, setup.N0, setup.T0)
                @test isapprox(c * est1.estimate, est2.estimate, rtol=1e-4)
                # Weights should be the same
                @test isapprox(est1.weights.omega, est2.weights.omega, atol=1e-4)
                @test isapprox(est1.weights.lambda, est2.weights.lambda, atol=1e-4)
            end
        end
    end

    @testset "Block shift tests" begin
        rng = MersenneTwister(42)
        setup = random_low_rank(rng=rng)
        Y = setup.Y
        N0 = setup.N0
        T0 = setup.T0
        N = size(Y, 1)
        T = size(Y, 2)
        exposed = (1:N) .> N0

        for c in [1e-6, 0.25, 1e6]
            # Shift exposed post: increases τ by c
            Y1 = copy(Y)
            Y1[exposed, (T0+1):T] .+= c
            for estimator in [synthdid_estimate, sc_estimate, did_estimate]
                est0 = estimator(Y, N0, T0)
                est1 = estimator(Y1, N0, T0)
                @test isapprox(est1.estimate, est0.estimate + c, rtol=1e-6)
            end

            # Shift exposed pre: decreases τ by c (not SC)
            Y2 = copy(Y)
            Y2[exposed, 1:T0] .+= c
            for estimator in [synthdid_estimate, did_estimate]
                est0 = estimator(Y, N0, T0)
                est2 = estimator(Y2, N0, T0)
                @test isapprox(est2.estimate, est0.estimate - c, atol=1e-10 + 1e-6 * abs(est0.estimate))
            end

            # Shift control pre: increases τ by c (not SC)
            Y3 = copy(Y)
            Y3[.!exposed, 1:T0] .+= c
            for estimator in [synthdid_estimate, did_estimate]
                est0 = estimator(Y, N0, T0)
                est3 = estimator(Y3, N0, T0)
                @test isapprox(est3.estimate, est0.estimate + c, atol=1e-10 + 1e-6 * abs(est0.estimate))
            end

            # Shift control post: decreases τ by c (not SC)
            Y4 = copy(Y)
            Y4[.!exposed, (T0+1):T] .+= c
            for estimator in [synthdid_estimate, did_estimate]
                est0 = estimator(Y, N0, T0)
                est4 = estimator(Y4, N0, T0)
                @test isapprox(est4.estimate, est0.estimate - c, atol=1e-10 + 1e-6 * abs(est0.estimate))
            end
        end
    end

    @testset "Variance estimation" begin
        # California has N1=1, so only placebo SE works (jackknife/bootstrap need N1>1)
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ = synthdid_estimate(setup.Y, setup.N0, setup.T0)

        Random.seed!(123)
        v_plac = vcov(τ; method=:placebo, replications=20)
        @test v_plac > 0
        println("Placebo SE (California): ", round(sqrt(v_plac), digits=2))

        # Use random data with multiple treated units for jackknife/bootstrap
        rng = MersenneTwister(42)
        rsetup = random_low_rank(rng=rng)
        τ_r = synthdid_estimate(rsetup.Y, rsetup.N0, rsetup.T0)

        v_jk = vcov(τ_r; method=:jackknife)
        @test v_jk > 0
        @test sqrt(v_jk) < abs(τ_r.estimate)

        Random.seed!(123)
        v_boot = vcov(τ_r; method=:bootstrap, replications=20)
        @test v_boot > 0

        println("Jackknife SE (random):  ", round(sqrt(v_jk), digits=2))
        println("Bootstrap SE (random):  ", round(sqrt(v_boot), digits=2))
    end

    @testset "Show method" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ = synthdid_estimate(setup.Y, setup.N0, setup.T0)
        # Just verify it doesn't error
        buf = IOBuffer()
        show(buf, τ)
        output = String(take!(buf))
        @test contains(output, "Synthetic Diff-in-Diff")
        @test contains(output, "τ̂")
    end

    @testset "Panel matrices returns X and labels" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        @test haskey(setup, :X)
        @test size(setup.X, 3) == 0  # no covariates
        @test haskey(setup, :units)
        @test haskey(setup, :times)
        @test length(setup.units) == size(setup.Y, 1)
        @test length(setup.times) == size(setup.Y, 2)
        @test "California" in setup.units
    end

    @testset "Units and times in estimate" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ = synthdid_estimate(setup.Y, setup.N0, setup.T0;
                              units=setup.units, times=collect(setup.times))
        @test length(τ.setup.units) == size(setup.Y, 1)
        @test length(τ.setup.times) == size(setup.Y, 2)
    end

    @testset "synthdid_controls" begin
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ_sdid = synthdid_estimate(setup.Y, setup.N0, setup.T0;
                                   units=setup.units, times=collect(setup.times))
        τ_sc = sc_estimate(setup.Y, setup.N0, setup.T0;
                           units=setup.units, times=collect(setup.times))

        # Single estimate
        ctrl = synthdid_controls(τ_sdid; mass=0.9)
        @test hasproperty(ctrl, :unit)
        @test nrow(ctrl) >= 1
        @test nrow(ctrl) <= setup.N0

        # Multiple estimates
        ctrl2 = synthdid_controls([τ_sdid, τ_sc]; mass=0.95)
        @test size(ctrl2, 2) == 3  # :unit + 2 estimator columns

        # Time weights
        ctrl_t = synthdid_controls(τ_sdid; weight_type="lambda", mass=0.9)
        @test nrow(ctrl_t) >= 1
    end

    @testset "Covariate support (empty covariates)" begin
        # Passing X=nothing should give same result as without
        df = california_prop99()
        setup = panel_matrices(df, :State, :Year, :PacksPerCapita, :treated)
        τ1 = synthdid_estimate(setup.Y, setup.N0, setup.T0)
        τ2 = synthdid_estimate(setup.Y, setup.N0, setup.T0; X=nothing)
        @test isapprox(τ1.estimate, τ2.estimate, atol=1e-10)
        @test isempty(τ1.weights.beta)
        @test isempty(τ2.weights.beta)
    end
end
