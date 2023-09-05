mutable struct AR_Markov
    y::Matrix{Float64}
    p::Int
    ϕ::Vector{Float64}
    aic::Float64
    aicc::Float64
    var_resid::Float64
    resid::Vector{Float64}
    p_values::Vector{Float64}
    coef_table
    ols
    fitted_y::Vector{Float64}
    fitted_X::Matrix{Float64}
    function AR_Markov(y::Matrix{Float64}, p::Int)
        assert_series_without_missing(y)
        return new(y,
            p,
            zeros(Float64, p),
            zero(Float64),
            zero(Float64),
            zero(Float64),
            zeros(Float64, 1),
            zeros(Float64, p),
            nothing,
            nothing,
            zeros(Float64, 1),
            zeros(Float64, 1, 1)
        )
    end
end

aic(ar::AR_Markov) = ar.aic
aicc(ar::AR_Markov) = ar.aicc
coef(ar::AR_Markov) = ar.ϕ
residuals(ar::AR_Markov) = ar.resid
residuals_variance(ar::AR_Markov) = ar.var_resid

mutable struct PARpMarkov
    y::Matrix{Float64}
    y_normalized::Matrix{Float64}
    clusters::Vector{Vector{Int64}}
    μ_stage::Vector{Float64}
    σ_stage::Vector{Float64}
    candidate_AR_stage::Vector{Vector{Vector{AR_Markov}}} # stage, cluster, p
    best_AR_stage::Vector{Vector{AR_Markov}} # stage, cluster
    n_stages::Int
    n_clusters::Int
    n_scenarios::Int
    p_lim::Int
    information_criteria::String
    series_with_zeros::Bool

    function PARpMarkov(y::Matrix{Float64}, clusters::Vector{Vector{Int64}}, p_lim::Int; information_criteria::String = "aic")
        assert_series_without_missing(y)
        series_with_zeros = series_with_only_zeros(y)
        n_stages = size(y, 1)
        n_clusters = maximum(maximum(clusters))
        n_scenarios = size(y, 2)
        y_normalized, μ_stage, σ_stage = normalize_series(y, n_stages)
        candidate_AR_stage = Vector{Vector{Vector{AR_Markov}}}(undef, 0)
        best_AR_stage = Vector{Vector{AR_Markov}}(undef, 0)
        return new(y, 
                y_normalized, 
                clusters, 
                μ_stage, 
                σ_stage, 
                candidate_AR_stage,
                best_AR_stage,
                n_stages, 
                n_clusters, 
                n_scenarios,
                p_lim,
                information_criteria,
                series_with_zeros
            )
    end
end

num_scenarios(par::PARpMarkov) = par.n_scenarios
num_stages(par::PARpMarkov) = par.n_stages
num_clusters(par::PARpMarkov) = par.n_clusters
p_limit(par::PARpMarkov) = par.p_lim
get_cluster_scenarios(par::PARpMarkov, stage::Int, cluster::Int) = (par.clusters[stage] .== cluster)

function residuals_of_best_model_at_stage(par::PARpMarkov, stage::Int)
    return residuals(par.best_AR_stage[stage])
end
function residuals_of_best_models_at_stage(par_models::Vector{PARpMarkov}, current_stage_to_predict::Int)
    residuals = Vector{Vector{Float64}}(undef, 0)
    for pm in par_models
        push!(residuals, residuals_of_best_model_at_stage(pm, current_stage_to_predict))
    end
    return concatenate_from_the_bottom_elements(residuals)
end

function build_y_X(y_normalized::Matrix{Float64}, p::Int, stage::Int)
    y_normalized_at_correct_stage = y_normalized[stage, :]
    n_scenarios = length(y_normalized_at_correct_stage)
    X = Matrix{Float64}(undef, n_scenarios, p)
    for i in 1:p
        lag_stage = stage - i
        if lag_stage < 1
            X[:, i] .= NaN
        else
            X[:, i] .= y_normalized[lag_stage, :]
        end
    end
    return y_normalized_at_correct_stage, X
end

function fit_ar!(ar::AR_Markov; stage::Int = 1)
    y_to_fit, X_to_fit = build_y_X(ar.y, ar.p, stage)
    ar.fitted_X = X_to_fit
    ar.fitted_y = y_to_fit
    if all(iszero, y_to_fit)
        ar.ols = nothing
        ar.ϕ = [1.0]
        ar.coef_table = nothing
        ar.p_values = zeros(1)
        ar.resid = y_to_fit
        ar.var_resid = 0.0
        ar.aic = 0.0
        ar.aicc = 0.0
    else
        ols = lm(X_to_fit, y_to_fit)
        ar.ols = ols
        ar.ϕ = coef(ols)
        ar.coef_table = coeftable(ols)
        ar.p_values = ar.coef_table.cols[ar.coef_table.pvalcol]
        ar.resid = y_to_fit - X_to_fit * ar.ϕ
        ar.var_resid = var(ar.resid)
        n = length(y_to_fit) + ar.p
        # A different but similiar form of aic can be found on Akaike original paper.
        # A note on the difference between the general form -2L + 2k and this one 
        # can also be found on wikipeadia talking about least squares estimators
        ar.aic = n * log(ar.var_resid * (n - 1)) + 2 * ar.p
        ar.aicc = ar.aic + (2 * ar.p^2 + 2 * ar.p)/(n - ar.p - 1)
    end
    return ar
end

function fit_par!(par::PARpMarkov)
    # fit all AR models
    for stage in 1:num_stages(par)
        candiate_ar_per_stage = Vector{AR_Markov}[]
        for cluster in 1:num_clusters(par)
            scenarios = get_cluster_scenarios(par, stage, cluster)
            candiate_ar_per_stage_cluster = AR_Markov[]
            for p in 1:par.p_lim
                candidate_ar_at_stage_cluster = AR_Markov(par.y_normalized[:, scenarios], p)
                fit_ar!(candidate_ar_at_stage_cluster; stage = stage)
                push!(candiate_ar_per_stage_cluster, candidate_ar_at_stage_cluster)
            end
            push!(candiate_ar_per_stage, candiate_ar_per_stage_cluster)
        end
        push!(par.candidate_AR_stage, candiate_ar_per_stage)
        best_model = select_best_model.(candiate_ar_per_stage, par.information_criteria)
        push!(par.best_AR_stage, best_model)
    end
    return par
end

function assert_same_number_of_stages(par_models::Vector{PARpMarkov})
    @assert length(unique(num_stages.(par_models))) == 1
end
function assert_same_p_limit(par_models::Vector{PARpMarkov})
    @assert length(unique(p_limit.(par_models))) == 1
end

function build_ar_models(results::ClusteringResults; header::Vector{String}, num_stages::Int, num_clus::Int, path_case::String, p_lim::Int = 6, verbose::Bool = false)
    params_ar_per_stage = DataStructures.OrderedDict()
    for agent_idx in eachindex(header)
        if verbose
            println("Gauging station $(header[agent_idx])...")
        end
        # Fit model
        par = PARpMarkov(results.input[agent_idx, :, :], results.cluster_results, p_lim)
        fit_par!(par)
        # Build dict
        AR_models_per_stage = Vector{DataStructures.OrderedDict}(undef, num_stages)
        params_ar_per_stage[header[agent_idx]] = DataStructures.OrderedDict(
            "Gauging station" => header[agent_idx],
            "AR_models_per_stage" => AR_models_per_stage
        )
        for stage_idx in 1:num_stages
            AR_models_per_cluster = Vector{DataStructures.OrderedDict}(undef, num_clus)
            AR_models_per_stage[stage_idx] = DataStructures.OrderedDict(
                "Stage" => stage_idx,
                "AR_models_per_cluster" => AR_models_per_cluster
            )
            for cluster_idx in 1:num_clus
                ar = par.best_AR_stage[stage_idx][cluster_idx]
                AR_models_per_cluster[cluster_idx] = DataStructures.OrderedDict(
                    "Cluster" => cluster_idx,
                    "Model coefficients" => ar.ϕ,
                )
            end
        end
    end
    write_par_params(path_case, params_ar_per_stage)
    return params_ar_per_stage
end
