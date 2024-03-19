module PAR

using LinearAlgebra
using StatsBase
using Statistics
using GLM
using PSRClassesClassicInterface
using Clustering
using Random
using DataStructures
using JSON

import GLM.coef

export fit_par!
export simulate_par
export PARp
export PARpA
export PARpMarkov

include("utils.jl")
include("PARp.jl")
include("PARpA.jl")

include("PARpMarkov/MarkovClustering.jl")
include("PARpMarkov/PARpMarkov.jl")
include("PARpMarkov/outputs.jl")

end # module
