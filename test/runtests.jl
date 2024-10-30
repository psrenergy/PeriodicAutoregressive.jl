using PeriodicAutoregressive

using Aqua
using Test

include("aqua.jl")
include("PARp.jl")
include("PARpA.jl")
include("utils.jl")

# data
include("data/funil_grande.jl")
include("data/batalha.jl")
include("data/camargos.jl")
include("data/inflow_with_stages_0_variance.jl")

function test_all()
    @testset "Aqua.jl" begin
        test_aqua()
    end

    @testset "utils.jl" begin
        test_utils()
    end

    @testset "PAR(p)" begin
        test_PARp()
    end

    @testset "PAR(p)-A" begin
        test_PARpA()
    end

    return nothing
end

test_all()