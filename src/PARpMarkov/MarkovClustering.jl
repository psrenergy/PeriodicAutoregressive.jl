mutable struct ClusteringOptions
    data::PSRClassesClassicInterface.PSRClassesData
    path_case::String
    file_name::String
    header::Vector{String}
    num_clus::Int64
    num_stages::Int64
end

mutable struct ClusteringResults
    input::Array{Float64, 3}
    cluster_results::Vector{Vector{Int64}}
    markov_transition_matrices::Vector{Matrix{Float64}}
end

function read_variable_to_cluster(options::ClusteringOptions)

    full_path = joinpath(options.path_case, options.file_name)

    if isfile(full_path * ".bin")
        reader = PSRI.open(
            PSRClassesClassicInterface.GrafReader{PSRClassesClassicInterface.PSRIOGrafResultBinary}, 
            full_path;
            header = options.header
        )
    else
        reader = PSRI.open(
            PSRClassesClassicInterface.GrafReader{PSRClassesClassicInterface.PSRIOGrafResult}, 
            full_path;
            header = options.header
        )
    end

    input = filter_variable_to_cluster(reader, options)
    variable_to_cluster = get_ena(input, reader, options)

    return variable_to_cluster, input
end

function filter_variable_to_cluster(reader::PSRClassesClassicInterface.GrafReader{T}, options::ClusteringOptions) where T
    # Aggregate by blocks and filter in the selected stages
    variable_to_cluster = zeros(Float64, length(options.header), options.num_stages, PSRI.max_scenarios(reader))

    for stage in 1:options.num_stages
        for scenario = 1:PSRI.max_scenarios(reader)
            for block = 1:PSRI.max_blocks_stage(reader, stage)
                PSRI.goto(reader, stage, scenario, block)
                duration = PSRI.block_duration(options.data, stage, block) / PSRI.stage_duration(options.data, stage)
                variable_to_cluster[:, stage, scenario] .+= duration * reader.data
            end
        end
    end
    return variable_to_cluster
end

function cluster_variable_by_scenario(options::ClusteringOptions)
    variable_to_cluster, input = read_variable_to_cluster(options)
    cluster_results = Vector{Vector{Int64}}(undef, options.num_stages)
    for stage in axes(variable_to_cluster, 1)
        kmeans_result = kmeans(variable_to_cluster[stage:stage, :], options.num_clus; rng = MersenneTwister(1))
        cluster_results[stage] = sort_clusters(kmeans_result)
    end
    return cluster_results, input
end

function sort_clusters(clus::KmeansResult{Array{Float64,2},Float64,Int64})
    # Have centers sorted in a vector
    vec_centers = vec(clus.centers)
    old_order = collect(1:length(vec_centers))
    new_order = sortperm(vec_centers)

    new_assignments = replace(assignments(clus), (old_order .=> new_order)...)
    return new_assignments
end

function clustering(data::PSRClassesClassicInterface.PSRClassesData; path_case::String, file_name::String, header::Vector{String}, num_clus::Int64, num_stages::Int64)
    options = ClusteringOptions(data, path_case, file_name, header, num_clus, num_stages)
    cluster_results, input = cluster_variable_by_scenario(options)
    num_scenarios = size(input, 3)
    markov_transition_matrices = Vector{Matrix{Float64}}(undef, num_stages - 1)
    for stage in 1:num_stages - 1
        markov_transition_matrices[stage] = markov_transition_matrix(cluster_results, stage, stage + 1)
    end
    write_cluster_assignemts(cluster_results, num_stages = num_stages, num_scenarios = num_scenarios, path_case = path_case)
    write_transition_matrix(markov_transition_matrices, num_stages = num_stages, num_clus = num_clus, path_case = path_case)
    return ClusteringResults(
        input,
        cluster_results,
        markov_transition_matrices
    )
end

function markov_transition_matrix(cluster_results::Vector{Vector{Int64}}, from_stage::Int64, to_stage::Int64)
    # TODO - rewrite this function to be more clear what we are doing here
    states_from = cluster_results[from_stage]
    states_to = cluster_results[to_stage]

    @assert length(states_from) == length(states_to)

    num_cluster = maximum(states_from)
    matrix_prob = zeros(Float64, num_cluster, num_cluster)

    for (clus_from, clus_to) in zip(states_from, states_to)
        matrix_prob[clus_from, clus_to] += 1
    end

    # Normalize count by the number of observations in each state
    matrix_prob = matrix_prob ./ sum(matrix_prob, dims = 1)

    return matrix_prob
end

# -------------------------------------------------
# ENA helper functions
# -------------------------------------------------

function get_ena(inflow::Array{Float64, 3}, reader::PSRClassesClassicInterface.GrafReader{T}, options::ClusteringOptions; system::String = "SUDESTE") where T
    data = options.data
    num_scenarios = PSRI.max_scenarios(reader)

    # Read data
    nhydro = PSRI.max_elements(data, "PSRHydroPlant")
    iphydrosys = PSRI.get_map(data, "PSRHydroPlant", "PSRSystem")
    iphydrosts = PSRI.get_map(data, "PSRHydroPlant", "PSRGaugingStation")
    iphydTurbTo = PSRI.get_map(data, "PSRHydroPlant", "PSRHydroPlant", relation_type = PSRI.PMD.RELATION_TURBINE_TO)
    sysName = PSRI.get_name(data, "PSRSystem")
    hydroFPmed = PSRI.mapped_vector(data, "PSRHydroPlant", "FPMed", Float64)
    hydroExisting = PSRI.mapped_vector(data, "PSRHydroPlant", "Existing", Int32)
    PSRI.update_vectors!(data)

    sys_idx = findfirst(x->x==system, sysName)

    hydroFPsys = zeros(nhydro)
    for i in 1:nhydro
        hydroFPsys[i] = (1 .- hydroExisting[i]) * downstream_ρ_system(i, sys_idx, iphydrosys, hydroFPmed, iphydTurbTo)
    end

    ena = zeros(nhydro, options.num_stages, num_scenarios)
    for h in 1:nhydro, stage in 1:options.num_stages, scen in 1:num_scenarios
        ena[h, stage, scen] = get_ENA(inflow[iphydrosts[h], stage, scen], hydroFPsys[h])
    end

    ena_agg = sum(ena, dims = 1)[1, :, :]
    
    return ena_agg
end

function downstream_ρ_system(h::Int, sys::Int, hyd2sys::Vector{Int32}, fprodt::Vector{Float64}, turbto::Vector{Int32})

    rho = 0.0
    nextplant = h
    added = BitSet()
    
    while nextplant != 0
        if nextplant in added
            break
        else
            push!(added, nextplant)
        end
        
        if hyd2sys[nextplant] == sys
            rho += fprodt[nextplant]
        end
        
        nextplant = turbto[nextplant]
    end

    return rho
end

function get_ENA(inflow::Float64, FP::Float64)
    day_len = 24 # [h]
    ena = max(inflow * FP * day_len, 0.0)
    return ena # [m3/s] * [MW/ m3/s] * [h] = [MWh]
end
