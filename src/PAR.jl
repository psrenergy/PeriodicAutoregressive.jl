module PAR

using LinearAlgebra
using StatsBase
using Statistics
using GLM

import GLM.coef

export fit_par!
export simulate_par
export PARp
export PARpA

include("utils.jl")
include("PARp.jl")
include("PARpA.jl")

end # module
