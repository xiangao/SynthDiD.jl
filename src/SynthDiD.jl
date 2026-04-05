module SynthDiD

using LinearAlgebra
using Statistics
using DataFrames
using CSV
using Random
using Distributions
using RecipesBase

export synthdid_estimate, sc_estimate, did_estimate,
       panel_matrices, effect_curve, placebo,
       vcov, se, california_prop99, synthdid_controls,
       SynthDiDEstimate, SynthDiDWeights, SynthDiDSetup, SynthDiDOpts

# ─── Structs ────────────────────────────────────────────────────────

struct SynthDiDSetup
    Y::Matrix{Float64}
    N0::Int
    T0::Int
    X::Array{Float64, 3}          # N×T×C covariate array (from Gemini)
    units::Vector{Any}             # unit labels (from OpenCode)
    times::Vector{Any}             # time labels (from OpenCode)
end

struct SynthDiDWeights
    omega::Vector{Float64}
    lambda::Vector{Float64}
    beta::Vector{Float64}          # covariate coefficients (from Gemini)
    omega_intercept::Bool
    lambda_intercept::Bool
    vals::Vector{Float64}          # optimization objective trace
end

struct SynthDiDOpts
    zeta_omega::Float64
    zeta_lambda::Float64
    omega_intercept::Bool
    lambda_intercept::Bool
    update_omega::Bool
    update_lambda::Bool
    min_decrease::Float64
    max_iter::Int
end

struct SynthDiDEstimate
    estimate::Float64
    estimator::Symbol           # :synthdid, :sc, :did
    weights::SynthDiDWeights
    setup::SynthDiDSetup
    opts::SynthDiDOpts
end

# Make it behave like a number
Base.Float64(e::SynthDiDEstimate) = e.estimate
Base.convert(::Type{Float64}, e::SynthDiDEstimate) = e.estimate
Base.isless(e::SynthDiDEstimate, x::Real) = isless(e.estimate, x)
Base.isless(x::Real, e::SynthDiDEstimate) = isless(x, e.estimate)
Base.:-(e::SynthDiDEstimate) = -e.estimate
Base.:-(a::SynthDiDEstimate, b::SynthDiDEstimate) = a.estimate - b.estimate
Base.:-(a::SynthDiDEstimate, b::Real) = a.estimate - b
Base.:-(a::Real, b::SynthDiDEstimate) = a - b.estimate
Base.:+(a::SynthDiDEstimate, b::Real) = a.estimate + b
Base.:*(a::Real, b::SynthDiDEstimate) = a * b.estimate

function Base.show(io::IO, est::SynthDiDEstimate)
    N0 = est.setup.N0
    T0 = est.setup.T0
    N = size(est.setup.Y, 1)
    T = size(est.setup.Y, 2)
    name = Dict(:synthdid => "Synthetic Diff-in-Diff",
                :sc => "Synthetic Control",
                :did => "Diff-in-Diff")[est.estimator]
    println(io, "$name Estimate")
    println(io, "─" ^ 35)
    println(io, "  τ̂  = ", round(est.estimate, digits=4))
    if !isempty(est.weights.beta)
        println(io, "  β = ", round.(est.weights.beta, digits=4))
    end
    println(io, "  N₀ = $N0,  N₁ = $(N - N0)")
    println(io, "  T₀ = $T0,  T₁ = $(T - T0)")
    ω = est.weights.omega
    λ = est.weights.lambda
    N0_eff = sum(ω) > 0 ? 1 / sum(ω .^ 2) : N0
    T0_eff = sum(λ) > 0 ? 1 / sum(λ .^ 2) : T0
    println(io, "  N₀_eff = ", round(N0_eff, digits=1),
            ",  T₀_eff = ", round(T0_eff, digits=1))
end

# ─── Data ──────────────────────────────────────────────────────────

"""
    california_prop99()

Load the California Proposition 99 dataset (Abadie, Diamond & Hainmueller 2010).
Returns a DataFrame with columns: State, Year, PacksPerCapita, treated.
"""
function california_prop99()
    datapath = joinpath(@__DIR__, "..", "data", "california_prop99.csv")
    CSV.read(datapath, DataFrame; delim=';')
end

# ─── Panel utilities ───────────────────────────────────────────────

"""
    panel_matrices(panel, unit, time, outcome, treatment; covariates=Symbol[], treated_last=true)

Convert a long balanced panel DataFrame to matrix format.
Returns a NamedTuple (Y, N0, T0, W, X, units, times) where Y is N×T and X is N×T×C.
"""
function panel_matrices(panel::DataFrame, unit::Symbol, time::Symbol,
                        outcome::Symbol, treatment::Symbol;
                        covariates::Vector{Symbol}=Symbol[],
                        treated_last::Bool=true)
    all_cols = vcat([unit, time, outcome, treatment], covariates)
    df = select(panel, all_cols)
    @assert !any(ismissing, eachcol(df) |> Iterators.flatten) "Missing values in panel"

    # Check balanced panel
    units = sort(unique(df[!, unit]))
    times = sort(unique(df[!, time]))
    @assert nrow(df) == length(units) * length(times) "Panel must be balanced"

    # Check treatment is 0/1
    @assert all(x -> x in (0, 1), df[!, treatment]) "Treatment must be 0 or 1"

    # Sort by unit then time
    sort!(df, [unit, time])

    N = length(units)
    T = length(times)
    C = length(covariates)

    Y = Matrix{Float64}(reshape(df[!, outcome], T, N)')
    W = Matrix{Float64}(reshape(df[!, treatment], T, N)')
    X = zeros(N, T, C)
    for c in 1:C
        X[:, :, c] = Matrix{Float64}(reshape(df[!, covariates[c]], T, N)')
    end

    # Identify treated units and treatment start
    w = vec(any(W .== 1, dims=2))
    T0_idx = findfirst(vec(any(W .== 1, dims=1)))
    @assert !isnothing(T0_idx) "No treated units found"
    T0 = T0_idx - 1
    N0 = sum(.!w)

    # Check simultaneous adoption
    @assert all(W[.!w, :] .== 0) "Control units must never be treated"
    @assert all(W[:, 1:T0] .== 0) "No units treated before T0"
    @assert all(W[w, (T0+1):end] .== 1) "All treated units must be treated after T0"

    if treated_last
        order = sortperm(W[:, T0+1])
        Y = Y[order, :]
        W = W[order, :]
        for c in 1:C
            X[:, :, c] = X[order, :, c]
        end
    end

    # Store unit/time labels
    (Y=Y, N0=N0, T0=T0, W=W, X=X,
     units=treated_last ? units[sortperm(W[invperm(order), T0+1])] : units,
     times=times)
end

# ─── Core algorithm ────────────────────────────────────────────────

"""
    contract3(X, beta)

Compute the sum over the third dimension of X weighted by beta: ∑ X[:,:,c] * beta[c].
"""
function contract3(X::Array{Float64, 3}, beta::Vector{Float64})
    N, T, C = size(X)
    out = zeros(N, T)
    if isempty(beta) || C == 0
        return out
    end
    for c in 1:C
        out .+= beta[c] .* X[:, :, c]
    end
    return out
end

"""
Collapse Y to (N0+1) × (T0+1) by averaging treated rows and post-treatment columns.
"""
function collapsed_form(Y::Matrix{Float64}, N0::Int, T0::Int)
    N, T = size(Y)
    N1 = N - N0
    T1 = T - T0

    Yc = zeros(N0 + 1, T0 + 1)
    Yc[1:N0, 1:T0] = Y[1:N0, 1:T0]
    Yc[1:N0, T0+1] = vec(mean(Y[1:N0, (T0+1):T], dims=2))
    Yc[N0+1, 1:T0] = vec(mean(Y[(N0+1):N, 1:T0], dims=1))
    Yc[N0+1, T0+1] = mean(Y[(N0+1):N, (T0+1):T])
    return Yc
end

"""
    sc_weight_fw_covariates(Y, X, zeta_lambda, zeta_omega; kwargs...)

Alternating optimization for omega, lambda, and beta (covariate coefficients).
"""
function sc_weight_fw_covariates(Y::Matrix{Float64}, X::Array{Float64, 3},
                                 zeta_lambda::Float64, zeta_omega::Float64;
                                 lambda_intercept::Bool=true, omega_intercept::Bool=true,
                                 min_decrease::Float64=1e-3, max_iter::Int=1000,
                                 lambda::Union{Nothing,Vector{Float64}}=nothing,
                                 omega::Union{Nothing,Vector{Float64}}=nothing,
                                 beta::Union{Nothing,Vector{Float64}}=nothing,
                                 update_lambda::Bool=true, update_omega::Bool=true)
    N0_plus_1, T0_plus_1 = size(Y)
    N0 = N0_plus_1 - 1
    T0 = T0_plus_1 - 1
    C = size(X, 3)

    if isnothing(lambda); lambda = fill(1.0 / T0, T0); end
    if isnothing(omega);  omega  = fill(1.0 / N0, N0); end
    if isnothing(beta);   beta   = zeros(C); end

    function update_weights(Y_curr, λ_curr, ω_curr)
        # Update lambda
        Y_λ = lambda_intercept ? Y_curr[1:N0, :] .- mean(Y_curr[1:N0, :], dims=1) : Y_curr[1:N0, :]
        if update_lambda
            λ_curr = fw_step(Y_λ[:, 1:T0], λ_curr, Y_λ[:, T0+1], N0 * zeta_lambda^2)
        end
        err_λ = Y_λ * vcat(λ_curr, -1.0)

        # Update omega
        Y_ω = omega_intercept ? Y_curr[:, 1:T0]' .- mean(Y_curr[:, 1:T0]', dims=1) : Y_curr[:, 1:T0]'
        if update_omega
            ω_curr = fw_step(Y_ω[:, 1:N0], ω_curr, Y_ω[:, N0+1], T0 * zeta_omega^2)
        end
        err_ω = Y_ω * vcat(ω_curr, -1.0)

        val = zeta_omega^2 * dot(ω_curr, ω_curr) + zeta_lambda^2 * dot(λ_curr, λ_curr) +
              dot(err_ω, err_ω) / T0 + dot(err_λ, err_λ) / N0
        return (val=val, lambda=λ_curr, omega=ω_curr, err_lambda=err_λ, err_omega=err_ω)
    end

    vals = fill(NaN, max_iter)
    t = 0
    Y_beta = Y .- (C > 0 ? collapsed_form(contract3(X, beta), N0, T0) : 0.0)
    w = update_weights(Y_beta, lambda, omega)

    while t < max_iter
        t += 1
        if C > 0
            grad_beta = zeros(C)
            for c in 1:C
                Xc = collapsed_form(X[:, :, c], N0, T0)
                term1 = dot(w.err_lambda, Xc[1:N0, :] * vcat(w.lambda, -1.0)) / N0
                term2 = dot(w.err_omega, Xc[:, 1:T0]' * vcat(w.omega, -1.0)) / T0
                grad_beta[c] = -(term1 + term2)
            end
            alpha = 1.0 / t
            beta .-= alpha .* grad_beta
            Y_beta = Y .- collapsed_form(contract3(X, beta), N0, T0)
            w = update_weights(Y_beta, w.lambda, w.omega)
        end
        vals[t] = w.val
        if t >= 2 && abs(vals[t-1] - vals[t]) <= min_decrease^2
            break
        end
    end

    return (lambda=w.lambda, omega=w.omega, beta=beta, vals=vals)
end

"""
Default sparsification: zero out entries ≤ max/4, renormalize.
"""
function sparsify(v::Vector{Float64})
    threshold = maximum(v) / 4
    w = copy(v)
    w[w .<= threshold] .= 0.0
    s = sum(w)
    s > 0 ? w ./ s : w
end

"""
    fw_step(A, x, b, eta; alpha=nothing)

One Frank-Wolfe step for min_x ||Ax - b||² + η||x||² over the unit simplex.
If alpha is given, use fixed step size; otherwise use exact line search.
"""
function fw_step(A::AbstractMatrix{Float64}, x::Vector{Float64},
                 b::AbstractVector{Float64}, eta::Float64;
                 alpha::Union{Nothing,Float64}=nothing)
    Ax = A * x
    half_grad = A' * (Ax - b) + eta * x
    i = argmin(half_grad)

    if !isnothing(alpha)
        x_new = x .* (1 - alpha)
        x_new[i] += alpha
        return x_new
    else
        d_x = -copy(x)
        d_x[i] += 1.0
        if all(d_x .== 0)
            return copy(x)
        end
        d_err = A[:, i] - Ax
        step = -dot(half_grad, d_x) / (dot(d_err, d_err) + eta * dot(d_x, d_x))
        constrained_step = clamp(step, 0.0, 1.0)
        return x + constrained_step * d_x
    end
end

"""
    sc_weight_fw(Y, zeta; intercept=true, lambda=nothing, min_decrease=1e-3, max_iter=1000)

Frank-Wolfe solver for synthetic control weights with exact line search.
Y is N0 × (T0+1), where columns 1:T0 are donors and column T0+1 is the target.
"""
function sc_weight_fw(Y::AbstractMatrix{Float64}, zeta::Float64;
                      intercept::Bool=true,
                      lambda::Union{Nothing,Vector{Float64}}=nothing,
                      min_decrease::Float64=1e-3,
                      max_iter::Int=1000)
    T0 = size(Y, 2) - 1
    N0 = size(Y, 1)

    if isnothing(lambda)
        lambda = fill(1.0 / T0, T0)
    else
        lambda = copy(lambda)
    end

    Yw = if intercept
        Y .- mean(Y, dims=1)
    else
        copy(Y)
    end

    A = Yw[:, 1:T0]
    b = Yw[:, T0+1]
    eta = N0 * zeta^2

    vals = fill(NaN, max_iter)
    t = 0
    while t < max_iter
        t += 1
        lambda = fw_step(A, lambda, b, eta)
        err = Yw * vcat(lambda, -1.0)
        vals[t] = zeta^2 * dot(lambda, lambda) + dot(err, err) / N0
        if t >= 2 && vals[t-1] - vals[t] <= min_decrease^2
            break
        end
    end

    return (lambda=lambda, vals=vals)
end

# ─── Main estimators ────────────────────────────────────────────────

"""
    synthdid_estimate(Y, N0, T0; kwargs...)

Compute the synthetic diff-in-diff estimate (Algorithm 1 of Arkhangelsky et al. 2021).

# Arguments
- `Y`: N×T outcome matrix. Rows 1:N0 are control units, rows N0+1:N are treated.
- `N0`: number of control units.
- `T0`: number of pre-treatment periods.
- `X`: optional N×T×C covariate array.
- `eta_omega`: regularization parameter for omega. Default: (N1*T1)^(1/4).
- `eta_lambda`: regularization for lambda. Default: 1e-6.
- `omega_intercept`, `lambda_intercept`: whether to demean when estimating weights.
- `omega`, `lambda`, `beta`: optional pre-specified values.
- `min_decrease`: convergence tolerance.
- `max_iter`: maximum iterations.
- `do_sparsify`: whether to sparsify weights. Default: true.
- `max_iter_pre_sparsify`: iterations before sparsification. Default: 100.
- `units`, `times`: optional vectors of unit/time labels (used for plotting/diagnostics).
"""
function synthdid_estimate(Y::Matrix{Float64}, N0::Int, T0::Int;
                           X::Union{Nothing,Array{Float64, 3}} = nothing,
                           eta_omega::Float64 = ((size(Y,1) - N0) * (size(Y,2) - T0))^0.25,
                           eta_lambda::Float64 = 1e-6,
                           omega_intercept::Bool = true,
                           lambda_intercept::Bool = true,
                           omega::Union{Nothing,Vector{Float64}} = nothing,
                           lambda::Union{Nothing,Vector{Float64}} = nothing,
                           beta::Union{Nothing,Vector{Float64}} = nothing,
                           min_decrease::Float64 = 1e-5,
                           max_iter::Int = 10000,
                           do_sparsify::Bool = true,
                           max_iter_pre_sparsify::Int = 100,
                           zeta_omega::Union{Nothing,Float64} = nothing,
                           zeta_lambda::Union{Nothing,Float64} = nothing,
                           units = nothing,
                           times = nothing)
    N = size(Y, 1)
    T = size(Y, 2)
    @assert N > N0 && T > T0
    N1 = N - N0
    T1 = T - T0

    if isnothing(X); X = zeros(N, T, 0); end
    C = size(X, 3)
    if isnothing(beta); beta = zeros(C); end

    noise_level = std([diff(Y[i, 1:T0]) for i in 1:N0] |> x -> vcat(x...))
    if isnothing(zeta_omega);  zeta_omega = eta_omega * noise_level;  end
    if isnothing(zeta_lambda); zeta_lambda = eta_lambda * noise_level; end
    min_dec = isnothing(min_decrease) ? 1e-5 * noise_level : min_decrease

    update_omega = isnothing(omega)
    update_lambda = isnothing(lambda)

    if isnothing(omega); omega = fill(1.0 / N0, N0); end
    if isnothing(lambda); lambda = fill(1.0 / T0, T0); end

    vals = Float64[]

    if C == 0
        # ── No covariates: standard two-step FW ──
        # Estimate lambda (time weights)
        if update_lambda
            Yc = collapsed_form(Y, N0, T0)
            lambda_opt = sc_weight_fw(Yc[1:N0, :], zeta_lambda;
                                      intercept=lambda_intercept, lambda=lambda,
                                      min_decrease=min_dec,
                                      max_iter=do_sparsify ? max_iter_pre_sparsify : max_iter)
            if do_sparsify
                lambda_opt = sc_weight_fw(Yc[1:N0, :], zeta_lambda;
                                          intercept=lambda_intercept,
                                          lambda=sparsify(lambda_opt.lambda),
                                          min_decrease=min_dec, max_iter=max_iter)
            end
            lambda = lambda_opt.lambda
            vals = lambda_opt.vals
        end

        # Estimate omega (unit weights)
        if update_omega
            Yc = collapsed_form(Y, N0, T0)
            omega_opt = sc_weight_fw(Yc[:, 1:T0]', zeta_omega;
                                     intercept=omega_intercept, lambda=omega,
                                     min_decrease=min_dec,
                                     max_iter=do_sparsify ? max_iter_pre_sparsify : max_iter)
            if do_sparsify
                omega_opt = sc_weight_fw(Yc[:, 1:T0]', zeta_omega;
                                         intercept=omega_intercept,
                                         lambda=sparsify(omega_opt.lambda),
                                         min_decrease=min_dec, max_iter=max_iter)
            end
            omega = omega_opt.lambda
            if isempty(vals)
                vals = omega_opt.vals
            else
                vals = pairwise_sum_decreasing(vals, omega_opt.vals)
            end
        end
    else
        # ── With covariates: alternating optimization for ω, λ, β ──
        Yc = collapsed_form(Y, N0, T0)
        opt = sc_weight_fw_covariates(Yc, X, zeta_lambda, zeta_omega;
                                      lambda_intercept=lambda_intercept,
                                      omega_intercept=omega_intercept,
                                      min_decrease=min_dec, max_iter=max_iter,
                                      lambda=lambda, omega=omega, beta=beta,
                                      update_lambda=update_lambda, update_omega=update_omega)
        lambda = opt.lambda
        omega = opt.omega
        beta = opt.beta
        vals = opt.vals
    end

    # Compute estimate: τ = (-ω, 1/N1)' (Y - Xβ) (-λ, 1/T1)
    unit_weights = vcat(-omega, fill(1.0 / N1, N1))
    time_weights = vcat(-lambda, fill(1.0 / T1, T1))
    estimate = dot(unit_weights, (Y .- contract3(X, beta)) * time_weights)

    # Default labels if not provided
    if isnothing(units)
        units = vcat(["control_$i" for i in 1:N0], ["treated_$i" for i in 1:N1])
    end
    if isnothing(times)
        times = collect(1:T)
    end

    weights = SynthDiDWeights(omega, lambda, beta, omega_intercept, lambda_intercept, vals)
    setup = SynthDiDSetup(Y, N0, T0, X, units, times)
    opts = SynthDiDOpts(zeta_omega, zeta_lambda, omega_intercept, lambda_intercept,
                        update_omega, update_lambda, min_dec, max_iter)

    return SynthDiDEstimate(estimate, :synthdid, weights, setup, opts)
end

"""
    sc_estimate(Y, N0, T0; kwargs...)

Synthetic control estimator. Uses zero time weights and no unit intercept.
"""
function sc_estimate(Y::Matrix{Float64}, N0::Int, T0::Int;
                     eta_omega::Float64 = 1e-6, kwargs...)
    synthdid_estimate(Y, N0, T0;
                      eta_omega=eta_omega,
                      lambda=zeros(T0),
                      omega_intercept=false,
                      kwargs...)
end

"""
    did_estimate(Y, N0, T0; kwargs...)

Standard diff-in-diff estimator. Uses uniform weights.
"""
function did_estimate(Y::Matrix{Float64}, N0::Int, T0::Int; kwargs...)
    synthdid_estimate(Y, N0, T0;
                      omega=fill(1.0 / N0, N0),
                      lambda=fill(1.0 / T0, T0),
                      kwargs...)
end

# ─── Effect curve & placebo ─────────────────────────────────────────

"""
    effect_curve(est::SynthDiDEstimate)

Return the per-period treatment effect curve that was averaged to produce the estimate.
"""
function effect_curve(est::SynthDiDEstimate)
    Y = est.setup.Y
    X = est.setup.X
    N0 = est.setup.N0
    T0 = est.setup.T0
    N1 = size(Y, 1) - N0
    T1 = size(Y, 2) - T0
    ω = est.weights.omega
    λ = est.weights.lambda
    β = est.weights.beta

    Y_adj = Y .- contract3(X, β)
    unit_weights = vcat(-ω, fill(1.0 / N1, N1))
    tau_sc = unit_weights' * Y_adj          # 1×T vector
    tau_curve = tau_sc[(T0+1):end] .- dot(tau_sc[1:T0], λ)
    return vec(tau_curve)
end

"""
    placebo(est::SynthDiDEstimate; treated_fraction=nothing)

Compute a placebo estimate using pre-treatment data only.
"""
function placebo(est::SynthDiDEstimate; treated_fraction::Union{Nothing,Float64}=nothing)
    Y = est.setup.Y
    X = est.setup.X
    N0 = est.setup.N0
    T0 = est.setup.T0
    opts = est.opts
    weights = est.weights

    if isnothing(treated_fraction)
        treated_fraction = 1 - T0 / size(Y, 2)
    end
    placebo_T0 = floor(Int, T0 * (1 - treated_fraction))

    estimator = Dict(:synthdid => synthdid_estimate,
                     :sc => sc_estimate,
                     :did => did_estimate)[est.estimator]

    return estimator(Y[:, 1:T0], N0, placebo_T0;
                     X = size(X, 3) == 0 ? nothing : X[:, 1:T0, :],
                     beta = isempty(weights.beta) ? nothing : weights.beta,
                     zeta_omega = opts.zeta_omega,
                     zeta_lambda = opts.zeta_lambda,
                     omega_intercept = opts.omega_intercept,
                     lambda_intercept = opts.lambda_intercept,
                     min_decrease = opts.min_decrease,
                     max_iter = opts.max_iter)
end

# ─── Variance estimation ───────────────────────────────────────────

"""
    vcov(est::SynthDiDEstimate; method=:bootstrap, replications=200)

Estimate variance. Methods: :bootstrap, :jackknife, :placebo.
"""
function vcov(est::SynthDiDEstimate;
              method::Symbol=:bootstrap, replications::Int=200)
    if method == :bootstrap
        se_val = bootstrap_se(est, replications)
    elseif method == :jackknife
        se_val = jackknife_se(est)
    elseif method == :placebo
        se_val = placebo_se(est, replications)
    else
        error("Unknown method: $method. Use :bootstrap, :jackknife, or :placebo.")
    end
    return se_val^2
end

"""
    se(est::SynthDiDEstimate; method=:bootstrap, replications=200)

Standard error of the estimate.
"""
function se(est::SynthDiDEstimate; kwargs...)
    sqrt(vcov(est; kwargs...))
end

function sum_normalize(x::Vector{Float64})
    s = sum(x)
    s != 0 ? x ./ s : fill(1.0 / length(x), length(x))
end

# Bootstrap SE: Algorithm 2
function bootstrap_se(est::SynthDiDEstimate, replications::Int)
    samples = bootstrap_sample(est, replications)
    sqrt((replications - 1) / replications) * std(samples)
end

function bootstrap_sample(est::SynthDiDEstimate, replications::Int)
    Y = est.setup.Y
    X = est.setup.X
    N0 = est.setup.N0
    T0 = est.setup.T0
    N = size(Y, 1)
    opts = est.opts
    ω = est.weights.omega
    β = est.weights.beta

    if N0 == N - 1
        return fill(NaN, replications)
    end

    estimates = Vector{Float64}(undef, replications)
    count = 0
    while count < replications
        ind = sort(rand(1:N, N))
        n0_boot = sum(ind .<= N0)
        if n0_boot == 0 || n0_boot == N
            continue
        end
        omega_boot = sum_normalize(ω[sort(ind[ind .<= N0])])
        try
            est_boot = synthdid_estimate(Y[ind, :], n0_boot, T0;
                                         X = size(X, 3) == 0 ? nothing : X[ind, :, :],
                                         omega=omega_boot,
                                         beta = isempty(β) ? nothing : β,
                                         zeta_omega=opts.zeta_omega,
                                         zeta_lambda=opts.zeta_lambda,
                                         omega_intercept=opts.omega_intercept,
                                         lambda_intercept=opts.lambda_intercept,
                                         min_decrease=opts.min_decrease,
                                         max_iter=opts.max_iter)
            count += 1
            estimates[count] = est_boot.estimate
        catch
            continue
        end
    end
    return estimates
end

# Jackknife SE: Algorithm 3 (fixed-weights)
function jackknife_se(est::SynthDiDEstimate)
    Y = est.setup.Y
    X = est.setup.X
    N0 = est.setup.N0
    T0 = est.setup.T0
    N = size(Y, 1)
    opts = est.opts
    ω = est.weights.omega
    β = est.weights.beta

    if N0 == N - 1 || sum(ω .!= 0) == 1
        return NaN
    end

    u = zeros(N)
    for i in 1:N
        ind = setdiff(1:N, i)
        n0_jk = sum(ind .<= N0)
        omega_jk = sum_normalize(ω[ind[ind .<= N0]])
        est_jk = synthdid_estimate(Y[ind, :], n0_jk, T0;
                                    X = size(X, 3) == 0 ? nothing : X[ind, :, :],
                                    omega=omega_jk,
                                    lambda=copy(est.weights.lambda),
                                    beta = isempty(β) ? nothing : β,
                                    zeta_omega=opts.zeta_omega,
                                    zeta_lambda=opts.zeta_lambda,
                                    omega_intercept=opts.omega_intercept,
                                    lambda_intercept=opts.lambda_intercept,
                                    min_decrease=opts.min_decrease,
                                    max_iter=opts.max_iter)
        u[i] = est_jk.estimate
    end

    sqrt(((N - 1) / N) * (N - 1) * var(u))
end

# Placebo SE: Algorithm 4
function placebo_se(est::SynthDiDEstimate, replications::Int)
    Y = est.setup.Y
    X = est.setup.X
    N0 = est.setup.N0
    T0 = est.setup.T0
    N1 = size(Y, 1) - N0
    opts = est.opts
    ω = est.weights.omega
    β = est.weights.beta

    if N0 <= N1
        error("Must have more controls than treated units for placebo SE")
    end

    estimates = zeros(replications)
    for r in 1:replications
        ind = sample(1:N0, N0, replace=false)
        n0_plac = length(ind) - N1
        omega_plac = sum_normalize(ω[ind[1:n0_plac]])
        est_plac = synthdid_estimate(Y[ind, :], n0_plac, T0;
                                      X = size(X, 3) == 0 ? nothing : X[ind, :, :],
                                      omega=omega_plac,
                                      beta = isempty(β) ? nothing : β,
                                      zeta_omega=opts.zeta_omega,
                                      zeta_lambda=opts.zeta_lambda,
                                      omega_intercept=opts.omega_intercept,
                                      lambda_intercept=opts.lambda_intercept,
                                      min_decrease=opts.min_decrease,
                                      max_iter=opts.max_iter)
        estimates[r] = est_plac.estimate
    end
    sqrt((replications - 1) / replications) * std(estimates)
end

# ─── Diagnostics ───────────────────────────────────────────────────

"""
    synthdid_controls(estimates; sort_by=1, mass=0.9, weight_type="omega")

Return a DataFrame of control weights, sorted by magnitude.
- `weight_type`: "omega" for unit weights, "lambda" for time weights
- `mass`: retains controls with cumulative weight ≥ mass (default 90%)
- `sort_by`: which estimate to sort by (1-indexed)
"""
function synthdid_controls(estimates::Vector{SynthDiDEstimate};
                           sort_by::Int=1,
                           mass::Float64=0.9,
                           weight_type::String="omega")
    if weight_type ∉ ("omega", "lambda")
        error("weight_type must be \"omega\" or \"lambda\"")
    end

    if weight_type == "omega"
        weights_matrix = hcat([est.weights.omega for est in estimates]...)
        names = estimates[1].setup.units[1:estimates[1].setup.N0]
    else
        weights_matrix = hcat([est.weights.lambda for est in estimates]...)
        names = estimates[1].setup.times[1:estimates[1].setup.T0]
    end

    sorted_idx = sortperm(weights_matrix[:, sort_by], rev=true)
    sorted_weights = weights_matrix[sorted_idx, :]
    sorted_names = names[sorted_idx]

    cumsum_weights = cumsum(sorted_weights[:, sort_by])
    keep_idx = findfirst(x -> x >= mass, cumsum_weights)
    if isnothing(keep_idx)
        keep_idx = length(cumsum_weights)
    end

    result = sorted_weights[1:keep_idx, :]
    est_names = Symbol.([est.estimator for est in estimates])

    df = DataFrame(result, est_names; makeunique=true)
    insertcols!(df, 1, :unit => sorted_names[1:keep_idx])

    return df
end

function synthdid_controls(est::SynthDiDEstimate; kwargs...)
    synthdid_controls([est]; kwargs...)
end

# ─── Helpers ───────────────────────────────────────────────────────

"""
Component-wise sum of decreasing vectors, where NaN means the vector stopped decreasing.
"""
function pairwise_sum_decreasing(x::Vector{Float64}, y::Vector{Float64})
    n = max(length(x), length(y))
    result = fill(NaN, n)

    # Extend shorter vector with its last non-NaN value
    x_ext = fill(NaN, n)
    y_ext = fill(NaN, n)
    x_ext[1:length(x)] = x
    y_ext[1:length(y)] = y

    last_x = NaN
    last_y = NaN
    for i in 1:n
        if !isnan(x_ext[i]); last_x = x_ext[i]; end
        if !isnan(y_ext[i]); last_y = y_ext[i]; end

        xi = isnan(x_ext[i]) ? last_x : x_ext[i]
        yi = isnan(y_ext[i]) ? last_y : y_ext[i]

        if isnan(xi) && isnan(yi)
            result[i] = NaN
        else
            result[i] = (isnan(xi) ? 0.0 : xi) + (isnan(yi) ? 0.0 : yi)
        end
    end
    return result
end

# ─── Convenience: pass keyword args through sc/did to synthdid ──────

# Allow zeta_omega/zeta_lambda as direct kwargs (for bootstrap reuse)
function synthdid_estimate(Y::Matrix{Float64}, N0::Int, T0::Int,
                           omega::Vector{Float64}, lambda::Vector{Float64};
                           kwargs...)
    synthdid_estimate(Y, N0, T0; omega=omega, lambda=lambda, kwargs...)
end

# ─── Plot Recipes ───────────────────────────────────────────────────

"""
    plot(est::SynthDiDEstimate; overlay=0.0)

Parallel trends plot: treated vs synthetic control trajectory,
with vertical treatment onset line and time weight bars.
- `overlay`: slide control trajectory toward treated (0-1) to check parallel trends.
"""
@recipe function f(est::SynthDiDEstimate;
                   overlay=0.0)
    Y = est.setup.Y
    N0 = est.setup.N0
    T0 = est.setup.T0
    N = size(Y, 1)
    T_total = size(Y, 2)
    ω = est.weights.omega
    λ = est.weights.lambda
    N1 = N - N0

    # Construct synthetic control
    synth = vec(ω' * Y[1:N0, :])
    treated = vec(mean(Y[(N0+1):N, :], dims=1))
    times = 1:T_total

    # Apply overlay: interpolate between synthetic and treated
    if overlay > 0
        synth = synth .+ overlay .* (treated .- synth)
    end

    name = Dict(:synthdid => "Synthetic Diff-in-Diff",
                :sc => "Synthetic Control",
                :did => "Diff-in-Diff")[est.estimator]

    layout --> (2, 1)
    size --> (700, 500)
    link --> :x

    # Panel 1: Trajectories
    @series begin
        subplot := 1
        title --> name
        ylabel --> "Outcome"
        label := "Synthetic Control"
        color := :steelblue
        linewidth := 2
        times, synth
    end

    @series begin
        subplot := 1
        label := "Treated"
        color := :firebrick
        linewidth := 2
        times, treated
    end

    # Treatment onset vertical line
    @series begin
        subplot := 1
        seriestype := :vline
        label := ""
        color := :gray40
        linestyle := :dash
        linewidth := 1
        [T0 + 0.5]
    end

    # Treatment effect arrow annotation (post-treatment gap)
    post_synth = mean(synth[(T0+1):end])
    post_treat = mean(treated[(T0+1):end])

    @series begin
        subplot := 1
        seriestype := :path
        label := "τ̂ = $(round(est.estimate, digits=1))"
        color := :black
        linewidth := 2
        arrow := true
        [T0 + (T_total - T0) / 2, T0 + (T_total - T0) / 2],
        [post_synth, post_treat]
    end

    # Panel 2: Time weights (bar chart)
    @series begin
        subplot := 2
        seriestype := :bar
        xlabel --> "Period"
        ylabel --> "Weight λ"
        label := ""
        color := :steelblue
        alpha := 0.7
        bar_width := 0.8
        times[1:T0], λ
    end

    # Treatment onset line in lower panel
    @series begin
        subplot := 2
        seriestype := :vline
        label := ""
        color := :gray40
        linestyle := :dash
        linewidth := 1
        [T0 + 0.5]
    end
end

"""
    plot(estimates::Vector{SynthDiDEstimate})

Side-by-side comparison of multiple estimators (like R's synthdid_plot with facets).
"""
@recipe function f(estimates::Vector{SynthDiDEstimate})
    n = length(estimates)
    Y = estimates[1].setup.Y
    N0 = estimates[1].setup.N0
    T0 = estimates[1].setup.T0
    N = size(Y, 1)
    T_total = size(Y, 2)
    N1 = N - N0
    treated = vec(mean(Y[(N0+1):N, :], dims=1))
    times = 1:T_total

    layout --> (1, n)
    size --> (350 * n, 350)
    link --> :y

    for (col, est) in enumerate(estimates)
        ω = est.weights.omega
        synth = vec(ω' * Y[1:N0, :])

        name = Dict(:synthdid => "SDiD",
                    :sc => "SC",
                    :did => "DiD")[est.estimator]

        @series begin
            subplot := col
            title --> "$name (τ̂=$(round(est.estimate, digits=1)))"
            ylabel --> (col == 1 ? "Outcome" : "")
            xlabel --> "Period"
            label := (col == 1 ? "Synthetic" : "")
            color := :steelblue
            linewidth := 2
            times, synth
        end

        @series begin
            subplot := col
            label := (col == 1 ? "Treated" : "")
            color := :firebrick
            linewidth := 2
            times, treated
        end

        @series begin
            subplot := col
            seriestype := :vline
            label := ""
            color := :gray40
            linestyle := :dash
            linewidth := 1
            [T0 + 0.5]
        end
    end
end

end # module
