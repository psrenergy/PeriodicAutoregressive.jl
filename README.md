# PeriodicAutoregressive.jl
Repository with Periodic Autoregressive models.

This repository provides an implementation of **Periodic Autoregressive (PAR)** models, along with simple variations, designed for time series analysis in periodic data. PAR models are valuable tools in generating synthetic energy and flow scenarios, particularly in energy operation planning.

The full methodology behind PAR models is detailed in the article by Maceira et al. (2006), _Geração de Cenários Sintéticos de Energia e Vazão para o Planejamento da Operação Energética_, published in *Cadernos do IME: Série Estatística*. Access the original article in Portuguese [here](https://www.e-publicacoes.uerj.br/index.php/cadest/article/download/15760/11931).

To add the package you can do:

```julia
add PeriodicAutoregressive
```

and to run an example you can do:

```julia
using PeriodicAutoregressive
funil_grande = include(joinpath(pkgdir(PeriodicAutoregressive), "test", "data", "funil_grande.jl"))
batalha = include(joinpath(pkgdir(PeriodicAutoregressive), "test", "data", "batalha.jl"))

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
