module PAR

using LinearAlgebra
using StatsBase
using Statistics

export fit_par!
export simulate_par
export PARp

include("utils.jl")
include("PARp.jl")

end # module
