function write_par_params(path_case::String, dict::DataStructures.OrderedDict, ident = 4) 
    file_path = joinpath(path_case, "par_params.json")
    open(file_path, "w") do f
        JSON.print(f, dict, ident)
    end
    return file_path
end

function write_graf(path_case::String, file_name::String, blocks::Int, scenarios::Int, stages::Int, agents::Vector{String}, data::Array{Float64, 4})
    full_path = joinpath(path_case, file_name)
    io = PSRI.open(
        PSRClassesClassicInterface.GrafWriter{PSRClassesClassicInterface.PSRIOGrafResult},
        full_path;
        blocks = blocks,
        scenarios = scenarios,
        stages = stages,
        agents = agents,
        unit = "",
        # initial_stage = 1,
        # initial_year = 2023,
    )
    for t in 1:stages, s in 1:scenarios, b in 1:blocks
        PSRI.write_registry(
            io,
            data[:, t, s, b],
            t,
            s,
            b,
        )
    end
    PSRI.close(io)
    return full_path*".csv"
end

function write_cluster_assignemts(cluster_results::Vector{Vector{Int64}}; num_stages::Int, num_scenarios::Int, path_case::String)
    data = Array{Float64, 4}(undef, 1, num_stages, num_scenarios, 1)
    for stage in 1:num_stages
        data[:, stage, :, :] = cluster_results[stage]
    end
    return write_graf(path_case, "cluster_map", 1, num_scenarios, num_stages, ["Cluster"], data)
end

function write_transition_matrix(markov_transition_matrices::Vector{Matrix{Float64}}; num_stages::Int, num_clus::Int, path_case::String)
    data = Array{Float64, 4}(undef, num_clus, num_stages - 1, 1, num_clus)
    for stage in 1:(num_stages - 1)
        for block in 1:num_clus
            data[:, stage, :, block] = markov_transition_matrices[stage][:, block]
        end
    end
    header = ["$x" for x in 1:num_clus]
    return write_graf(path_case, "transition_matrix", num_clus, 1, num_stages - 1, header, data)
end

function write_backwards(path_case::String, agents::Vector{String}, data::Array{Float64, 4})
    # agent, stage, scen_forw, scen_back
    num_stages = size(data, 2)
    num_scen_forw = size(data, 3)
    num_scen_back = size(data, 4)
    return write_graf(path_case, "back", num_scen_back, num_scen_forw, num_stages, agents, data)
end
