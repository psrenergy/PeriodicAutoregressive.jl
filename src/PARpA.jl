mutable struct AR_A
    y::Vector{Float64}
    y_anual::Vector{Float64}
    p::Int
    ϕ::Vector{Float64}
    ϕ_A::Float64
    aic::Float64
    aicc::Float64
    var_resid::Float64
    resid::Vector{Float64}
    p_values::Vector{Float64}
    coef_table
    ols
    function AR_A(y::Vector{Float64}, y_anual::Vector{Float64}, p::Int)
        assert_series_without_missing(y)
        return new(y,
            y_anual,
            p,
            zeros(Float64, p),
            zero(Float64),
            zero(Float64),
            zero(Float64),
            zero(Float64),
            zeros(Float64, 2),
            zeros(Float64, p),
            nothing,
            nothing
        )
    end
end

aic(ar_a::AR_A) = ar_a.aic
aicc(ar_a::AR_A) = ar_a.aicc
coef(ar_a::AR_A) = [ar_a.ϕ; ar_a.ϕ_A]
residuals(ar_a::AR_A) = ar_a.resid
residuals_variance(ar_a::AR_A) = ar_a.var_resid

function get_y_anual(y::Vector{Float64})
    y_anual = zeros(length(y))
    for i in 13:length(y)
        y_anual[i] = mean(y[i-12:i-1])
    end
    return y_anual
end

mutable struct PARpA
    y::Vector{Float64}
    y_normalized::Vector{Float64}
    μ_stage::Vector{Float64}
    σ_stage::Vector{Float64}
    y_anual::Vector{Float64}
    y_anual_normalized::Vector{Float64}
    μ_anual::Vector{Float64}
    σ_anual::Vector{Float64}
    candidate_AR_A_stage::Vector{Vector{AR_A}}
    best_AR_A_stage::Vector{AR_A}
    seasonal::Int
    p_lim::Int
    information_criteria::String
    last_stage::Int

    function PARpA(y::Vector{Float64}, seasonal::Int, p_lim::Int; information_criteria::String = "aic")
        @assert 0 < p_lim < 12
        assert_series_without_missing(y)
        y_normalized, μ_stage, σ_stage = normalize_series(y, seasonal)
        y_anual = get_y_anual(y)
        # normalização da média anual
        y_anual_normalized, μ_anual, σ_anual = PAR.normalize_series(y_anual[13:end], 12)
        y_anual_normalized = vcat(NaN.*ones(12), y_anual_normalized)
        candidate_AR_A_stage = Vector{Vector{AR_A}}(undef, 0)
        best_AR_A_stage = Vector{AR_A}(undef, 0)
        last_stage = mod1(length(y), seasonal)
        return new(y, 
                y_normalized, 
                μ_stage, 
                σ_stage, 
                y_anual,
                y_anual_normalized,
                μ_anual, 
                σ_anual,
                candidate_AR_A_stage,
                best_AR_A_stage,
                seasonal, 
                p_lim,
                information_criteria,
                last_stage
            )
    end
end

num_stages(par::PARpA) = par.seasonal
p_limit(par::PARpA) = par.p_lim

function residuals_of_best_model_at_stage(par::PARpA, stage::Int)
    return residuals(par.best_AR_A_stage[stage])
end
function residuals_of_best_models_at_stage(par_models::Vector{PARpA}, current_stage_to_predict::Int)
    residuals = Vector{Vector{Float64}}(undef, 0)
    for pm in par_models
        push!(residuals, residuals_of_best_model_at_stage(pm, current_stage_to_predict))
    end
    return concatenate_from_the_bottom_elements(residuals)
end

## NEXT
function build_y_X(y_normalized::Vector{Float64}, y_anual_normalized::Vector{Float64}, 
                   p::Int, stage::Int, seasonal::Int)
    n_total = length(y_normalized)
    correct_stage_idx = collect(stage:seasonal:n_total)
    y_normalized_at_correct_stage = y_normalized[correct_stage_idx]
    y_anual_normalized_at_correct_stage = y_anual_normalized[correct_stage_idx]
    n_correct_stage = length(y_normalized_at_correct_stage)
    X = Matrix{Float64}(undef, n_correct_stage, p)
    for i in 1:p
        lag_stage_idx = correct_stage_idx .- i
        for (t, idx) in enumerate(lag_stage_idx)
            if idx < 1
                X[t, i] = NaN
            else
                X[t, i] = y_normalized[idx]
            end
        end
    end
    X_reg = hcat(X[p+1:end, :], y_anual_normalized_at_correct_stage[p+1:end])
    y_reg = y_normalized_at_correct_stage[p+1:end]
    return y_reg, X_reg
end

function fit_ar!(ar_A::AR_A; stage::Int = 1, par_seasonal::Int = 12)
    y_to_fit, X_to_fit = build_y_X(ar_A.y, ar_A.y_anual, ar_A.p, stage, par_seasonal)
    ols = lm(X_to_fit, y_to_fit)
    ar_A.ols = ols
    params = coef(ols)
    ar_A.coef_table = coeftable(ols)
    ar_A.p_values = ar_A.coef_table.cols[ar_A.coef_table.pvalcol]
    ar_A.ϕ = params[1:end-1]
    ar_A.ϕ_A = params[end]
    ar_A.resid = y_to_fit - X_to_fit * params
    ar_A.var_resid = var(ar_A.resid)
    n = length(y_to_fit) + ar_A.p
    # A different but similiar form of aic can be found on Akaike original paper.
    # A note on the difference between the general form -2L + 2k and this one 
    # can also be found on wikipeadia talking about least squares estimators
    ar_A.aic = n * log(ar_A.var_resid * (n - 1)) + 2 * (ar_A.p + 1)
    ar_A.aicc = ar_A.aic + (2 * (ar_A.p + 1)^2 + 2 * (ar_A.p+1))/(n - (ar_A.p + 1) - 1)
    return ar_A
end

function fit_par!(par::PARpA)
    # fit all AR models
    for stage in 1:num_stages(par)
        candiate_ar_a_per_stage = PAR.AR_A[]
        for p in 1:par.p_lim
            candidate_ar_a_at_stage = AR_A(par.y_normalized, par.y_anual_normalized, p)
            fit_ar!(candidate_ar_a_at_stage; stage = stage, par_seasonal = par.seasonal)
            push!(candiate_ar_a_per_stage, candidate_ar_a_at_stage)
        end
        push!(par.candidate_AR_A_stage, candiate_ar_a_per_stage)
        best_model = select_best_model(candiate_ar_a_per_stage, par.information_criteria)
        push!(par.best_AR_A_stage, best_model)
    end
    return par
end

function assert_same_number_of_stages(par_models::Vector{PARpA})
    @assert length(unique(num_stages.(par_models))) == 1
end
function assert_same_p_limit(par_models::Vector{PARpA})
    @assert length(unique(p_limit.(par_models))) == 1
end
function calculate_ϕ(pm::PARpA, stage::Int)
    n_stages = 12
    ϕ = zeros(pm.best_AR_A_stage[stage].p)
    for i in 1:pm.best_AR_A_stage[stage].p
        ϕ[i] = pm.best_AR_A_stage[stage].ϕ[i] * (pm.σ_stage[stage]/pm.σ_stage[mod1(stage - i, n_stages)])
    end
    return ϕ
end

function calculate_ψ(pm::PARpA, stage::Int)
    n_stages = 12
    return pm.best_AR_A_stage[stage].ϕ_A * (pm.σ_stage[stage] / pm.σ_anual[mod1(stage - 1, n_stages)])
end

simulate_par(par::PARpA, stepds_ahead::Int, n_scenarios::Int) = simulate_par([par], stepds_ahead, n_scenarios)
function simulate_par(par_models::Vector{PARpA}, steps_ahead::Int, n_scenarios::Int)
    assert_same_number_of_stages(par_models)
    assert_same_p_limit(par_models)
    n_stages = num_stages(par_models[1])
    n_models = length(par_models)
    scenarios_normalized = zeros(steps_ahead + n_stages, n_models, n_scenarios)
    scenarios = zeros(steps_ahead + n_stages, n_models, n_scenarios)
    # Fill the first part of scenarios with historical data
    for (i, pm) in enumerate(par_models)
        scenarios_normalized[1:n_stages, i,  :] .= pm.y_normalized[end-n_stages+1:end]
        scenarios[1:n_stages, i,  :] .= pm.y[end-n_stages+1:end]
    end
    # Simulate on the standardized series
    for t in 1:steps_ahead
        current_stage_to_predict = mod1(par_models[1].last_stage + t, n_stages)
        t_scen_idx = t + n_stages
        residuals_matrix = residuals_of_best_models_at_stage(par_models, current_stage_to_predict)
        cor_matrix = cor(residuals_matrix)
        ruido_normal = randn(n_scenarios, n_models) 
        ruido_correlacionado = ruido_normal * cholesky(cor_matrix).L
        for (i, pm) in enumerate(par_models)
            current_model_p = pm.best_AR_A_stage[current_stage_to_predict].p
            for s in 1:n_scenarios
                # Evaluate the deterministic parts of the scenarios
                autorregressive_normalized = dot(
                    scenarios_normalized[t_scen_idx-1 : -1 : t_scen_idx-current_model_p, i, s],
                    pm.best_AR_A_stage[current_stage_to_predict].ϕ
                )
                mean_last_12_ena = mean(scenarios[t_scen_idx-n_stages:t_scen_idx-1, i, s])
                anual_normalized = pm.best_AR_A_stage[current_stage_to_predict].ϕ_A * 
                                   (mean_last_12_ena - pm.μ_anual[mod1(current_stage_to_predict - 1, n_stages)]) / 
                                   pm.σ_anual[mod1(current_stage_to_predict - 1, n_stages)]

                # Calculate the noise and the parameters of the viable 3 parameter log normal
                lower_bound_log_normal = - (pm.μ_stage[current_stage_to_predict] / pm.σ_stage[current_stage_to_predict]) - 
                                            autorregressive_normalized - anual_normalized
                λ = (pm.best_AR_A_stage[current_stage_to_predict].var_resid/lower_bound_log_normal^2) + 1
                μ_log_normal = 0.5 * log(pm.best_AR_A_stage[current_stage_to_predict].var_resid/ (λ * (λ - 1)))
                σ_log_normal = sqrt(log(λ))
                ruido = exp(ruido_correlacionado[s, i]*σ_log_normal + μ_log_normal) + lower_bound_log_normal
                # Simulate scenarios
                scenarios_normalized[t_scen_idx, i, s] = autorregressive_normalized + anual_normalized + ruido
                scenarios[t_scen_idx, i, s] = scenarios_normalized[t_scen_idx, i, s] * pm.σ_stage[current_stage_to_predict] + pm.μ_stage[current_stage_to_predict]
            end
        end
    end
    return scenarios[n_stages+1:end, :, :]
end
