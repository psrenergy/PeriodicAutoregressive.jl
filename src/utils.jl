function assert_series_without_missing(y::Vector{Float64})
    for element in y
        isnan(element) && return throw(AssertionError("y cannot have missing values"))
    end
    return true
end

function series_with_only_zeros(y::Vector{Float64})
    for element in y 
        if !iszero(element)
            return false
        end
    end
    return true
end

function μ_σ_per_month(y::Vector{Float64}, seasonal::Int)
    μ = Vector{Float64}(undef, seasonal)
    σ = Vector{Float64}(undef, seasonal)
    for i in 1:seasonal
        seasonal_obs = y[i:seasonal:end]
        valid_obs = seasonal_obs[findall(!isnan, seasonal_obs)]
        μ[i] = mean(valid_obs)
        σ[i] = std(valid_obs)
    end
    return μ, σ
end

function normalize_series(y::Vector{Float64},
                          seasonal::Int)
    y_normalized = copy(y)
    μ, σ = μ_σ_per_month(y, seasonal)
    for i in 1:length(y_normalized)
        # We add 1e-5 to avoid dividing by 0.
        y_normalized[i] = (y_normalized[i] - μ[mod1(i, seasonal)]) / (σ[mod1(i, seasonal)] + 1e-5)
    end
    return y_normalized, μ, σ
end

"""
    concatenate_from_the_bottom_elements(vov::Vector{Vector{T}}) where T

Takes a vector of vectors and arrange the elements in a matrix based on the bottom element.

i.e.

a vector [1,2,3,4] and a vector [1, 2] should become the matrix [3 1; 4 2]
"""
function concatenate_from_the_bottom_elements(vov::Vector{Vector{T}}) where T
    n_series = length(vov)
    min_length = minimum(length.(vov))
    mat = zeros(min_length, n_series)
    for col in 1:n_series
        mat[:, col] = vov[col][end-min_length+1:end]
    end
    return mat
end

"""
    select_best_model(candidate_models::Vector, information_criteria::String)

Apply the information criteria to choose the best model.
"""
function select_best_model(candidate_models::Vector, information_criteria::String)
    if information_criteria == "aic"
        candidate_aic = map(aic, candidate_models)
        _, best_model_idx = findmin(candidate_aic)
        return candidate_models[best_model_idx]
    elseif information_criteria == "aicc"
        candidate_aicc = map(aicc, candidate_models)
        _, best_model_idx = findmin(candidate_aicc)
        return candidate_models[best_model_idx]
    elseif information_criteria == "fixed_at_p_lim"
        # TODO 
        # When the information criteria is to choose a fixed p for everybody 
        # the best model is always the one of order p_lim.
        # This implementation is very naive and can be optimzed (bu not fitting all models
        # from 1 to p_lim)
        return candidate_models[end]
    end
    return error()
end