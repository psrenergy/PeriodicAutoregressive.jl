# This function does an AR in the scenarios 
# TODO better document the procedure
function fit_ar_per_stage(
    forward::Array{Float64, 3}, 
    agent::Int, 
    stage::Int, 
    agent_key::String, # TODO we shound not need key and agent index
    observations::DataStructures.OrderedDict
)
    num_scenarios = size(forward, 2)
    y = forwards[agent, :, stage]

    aic_candidates = zeros(order)
    phi_candidates = Vector{Vector{Float64}}(undef, order)
    yhat_candidates = Vector{Vector{Float64}}(undef, order)
    residuals_candidates = Vector{Vector{Float64}}(undef, order)

    for order in 1:6 # Could be read from a parameter
        n_obs_lags = num_observation_lags(order, stage)
        n_forw_lags = num_forward_lags(order, stage)

        X = Matrix{Float64}(undef, num_scenarios, order)
        # TODO add intercept
        X[:,1:n_forw_lags] = forwards[agent, :, (stage - n_forw_lags):(stage - 1)]
        if n_obs_lags > 0
            for i = 1:n_obs_lags
                X[:, num_forward_lags + i] .= observations[agent_key]["lag $(i)"]
            end
        end
        phi = X \ y
        yhat = X * phi
        residuals = yhat .- y
        var_resid = var(residuals)

        # A different but similiar form of aic can be found on Akaike original paper.
        # A note on the difference between the general form -2L + 2k and this one 
        # can also be found on wikipeadia talking about least squares estimators
        n = num_scenarios
        aic_candidates[order] = n * log(var_resid * (n - 1)) + 2 * order
        phi_candidates[order] = phi
        yhat_candidates[order] = yhat
        residuals_candidates[order] = residuals
    end

    idx = argmin(aic_candidates)
    @assert !isnan(aic_candidates[idx])

    return phi_candidates[idx], yhat_candidates[idx], residuals_candidates[idx]
end

function num_observation_lags(order::Int, stage::Int)::Int
    if order >= stage
        return order - stage + 1
    else
        return 0
    end
end

function num_forward_lags(order::Int, stage::Int)::Int
    if stage > order
        return order
    else
        return order - obs_lags
    end
end

function build_backward(rng, y_hat::Float64, scenario_forward::Int, residuals::Vector{Float64}, num_scen_back::Int64)
    all_possible_backwards = y_hat[scenario_forward] .+ residuals
    # pick n_back elements of all_possibilities
    backwards = rand(rng, all_possible_backwards, num_scen_back)
    # The first backward must be the same as the forward
    backwards[1] = y_hat[scenario_forward] + residuals[scenario_forward]
    return backwards
end

function estimate_one_ar_per_stage(lags_ar_json_file::String, forwards::Array{Float64, 3}, num_scen_back::Int64)
    agents_historical_observations = JSON.parsefile(lags_ar_json_file; dicttype = DataStructures.OrderedDict)
    num_stages = size(forwards, 3)
    num_scen_forw = size(forwards, 2)
    num_agents = size(forwards, 1)

    parametros_ar_per_stage = DataStructures.OrderedDict()

    backward_ar_per_stage = Array{Float64}(undef, num_agents, num_scen_back, num_scen_forw, num_stages)
    agentKeys = collect(keys(agents_historical_observations))

    for agent = 1:num_agents
        agentKey = agentKeys[agent]
        AR_model_per_stage = Vector{DataStructures.OrderedDict}(undef, num_stages)
        parametros_ar_per_stage[agentKey] = DataStructures.OrderedDict(
            "mu_stage" => zeros(num_stages),
            "sigma_stage" => ones(num_stages),
            "AR_model_per_stage" => AR_model_per_stage
        )
        for stage = 1:num_stages
            phi, y_hat, residuals = fit_ar_per_stage(forward, agent, stage, agent_key, observations)
            AR_model_per_stage[stage] = DataStructures.OrderedDict(
                "coefs" => phi,
                "var_resid" => 0.0 # We don't use this field
            )
            for scen_forw = 1:num_scen_forw
                rng = Random.MersenneTwister(scen_forw)
                backward_ar_per_stage[agent, :, scenario_forward, stage] = build_backward(rng, y_hat, scen_forw, residuals, num_scen_back)
            end
        end
    end
    return parametros_ar_per_stage, backward_ar_per_stage
end