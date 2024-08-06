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

function series_with_all_same_observation(y::Vector{Float64})
    for i in 2:length(y)
        if y[i] != y[1]
            return false
        end
    end
    return true
end

function remove_columns_with_all_zeros(mat::Matrix{Float64})
    matrix_with_no_zeros = deepcopy(mat)
    indexes_of_removed_columns = Int[]
    indexes_of_kept_columns = Int[]
    column_new_index_map = Int[]
    for j in 1:size(mat, 2)
        if all(mat[:, j] .== 0)
            push!(indexes_of_removed_columns, j)
            push!(column_new_index_map, -1)
        else
            push!(column_new_index_map, j - length(indexes_of_removed_columns))
            push!(indexes_of_kept_columns, j)
        end
    end
    return matrix_with_no_zeros[:, indexes_of_kept_columns], indexes_of_removed_columns, column_new_index_map
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