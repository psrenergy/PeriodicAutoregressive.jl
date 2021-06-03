include(".\\src\\PAR.jl")

n_stages = 12
p_lim = 6
par_a_1 = PAR.PARpA(funil_grande, n_stages, p_lim; information_criteria = "aic");
par_1

using Statistics
k = 2
mean(camargos[k:k+11])

mean(camargos[1:12])