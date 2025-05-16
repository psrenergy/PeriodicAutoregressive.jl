# using PeriodicAutoregressive

# using Aqua
# using Test

# include("test_aqua.jl")
# include("PARp.jl")
# include("PARpA.jl")
# include("utils.jl")

# # data
# include("data/funil_grande.jl")
# include("data/batalha.jl")
# include("data/camargos.jl")
# include("data/inflow_with_stages_0_variance.jl")

# function test_all()
#     @testset "Aqua.jl" begin
#         test_aqua()
#     end

#     @testset "utils.jl" begin
#         test_utils()
#     end

#     @testset "PAR(p)" begin
#         test_PARp()
#     end

#     @testset "PAR(p)-A" begin
#         test_PARpA()
#     end

#     return nothing
# end

# test_all()

using Test

function test_modules(dir::AbstractString)
    result = String[]
    for (root, dirs, files) in walkdir(dir)
        append!(result, filter!(f -> occursin(r"test_(.)+\.jl", f), joinpath.(root, files)))
    end
    return result
end

for file in test_modules(@__DIR__)
    @testset "$(basename(file))" begin
        include(file)
    end
end
