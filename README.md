# PeriodicAutoregressive.jl

[![CI](https://github.com/psrenergy/PeriodicAutoregressive.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/psrenergy/PeriodicAutoregressive.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/psrenergy/PeriodicAutoregressive.jl/graph/badge.svg?token=7tA9ajgsLf)](https://codecov.io/gh/psrenergy/PeriodicAutoregressive.jl)

## Introduction

This repository provides an implementation of **Periodic Autoregressive (PAR)** models, along with simple variations, designed for time series analysis in periodic data. PAR models are valuable tools in generating synthetic energy and flow scenarios, particularly in energy operation planning.

The full methodology behind PAR models is detailed in the article by Maceira et al. (2006), _Geração de Cenários Sintéticos de Energia e Vazão para o Planejamento da Operação Energética_, published in *Cadernos do IME: Série Estatística*. Access the original article in Portuguese [here](https://typeset.io/pdf/geracao-de-cenarios-sinteticos-de-energia-e-vazao-para-o-34mgqextw5.pdf).

## Getting Started

### Installation

```julia
julia> ] add PeriodicAutoregressive
```

### Example: PAR(p) Model

```julia
using PeriodicAutoregressive

funil_grande = include(joinpath(pkgdir(PeriodicAutoregressive), "test", "data", "funil_grande.jl"))
batalha = include(joinpath(pkgdir(PeriodicAutoregressive), "test", "data", "batalha.jl"))

stages = 12
p_lim = 6

par_1 = PARp(funil_grande, stages, p_lim; information_criteria = "aic");
par_2 = PARp(batalha, stages, p_lim; information_criteria = "aicc");

fit_par!(par_1);
fit_par!(par_2);

steps_ahead = 100
number_of_scenarios = 1000

scenarios = simulate_par([par_1; par_2], steps_ahead, number_of_scenarios)
```
