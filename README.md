# PAR.jl
Repository with Periodic Autorregressive models.

This is an immature repository and should not be considered production-ready. The main goal of the repository is to implement simple variations of Periodic Autorregressive models.

To add the package you can do:
```julia
pkg> add https://github.com/psrenergy/PAR.jl.git
```

and to run an example you can do:
```julia
using PAR
funil_grande = include(joinpath(pkgdir(PAR), "test", "data", "funil_grande.jl"))
batalha = include(joinpath(pkgdir(PAR), "test", "data", "batalha.jl"))

n_stages = 12
p_lim = 6
par_1 = PARp(funil_grande, n_stages, p_lim; information_criteria = "aic");
par_2 = PARp(batalha, n_stages, p_lim; information_criteria = "aicc");
fit_par!(par_1);
fit_par!(par_2);

steps_ahead = 100
n_scenarios = 1000
scen = simulate_par([par_1; par_2], steps_ahead, n_scenarios)
```