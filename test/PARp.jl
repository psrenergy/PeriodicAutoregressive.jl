@testset "PAR(p)" begin
    n_stages = 12
    p_lim = 6
    par_1 = PARp(funil_grande, n_stages, p_lim; information_criteria = "aic");
    par_2 = PARp(batalha, n_stages, p_lim; information_criteria = "aicc");
    fit_par!(par_1);
    fit_par!(par_2);

    # TODO get all the estima files in this repository
    # Some random tests based on an estima.csv file
    @test par_1.μ_stage[2] ≈ 286.75 atol = 1e-2
    @test par_1.σ_stage[2] ≈ 124.45 atol = 1e-2
    @test par_1.best_AR_stage[4].p == 4

    @test par_2.best_AR_stage[2].ϕ[1] .≈ 0.66114 atol = 1e-2
    @test par_2.best_AR_stage[3].ϕ[1] .≈ 0.45045 atol = 1e-2

    # TODO test something about the simulations
    scen = simulate_par(par_1, 10, 100)
    scen = simulate_par([par_1; par_2], 100, 1000)

    p_lim = 3
    par_1_fixed_p = PARp(funil_grande, n_stages, p_lim; information_criteria = "fixed_at_p_lim");
    fit_par!(par_1_fixed_p);
    @test par_1_fixed_p.best_AR_stage[1].p == 3
    @test par_1_fixed_p.best_AR_stage[2].p == 3
    @test par_1_fixed_p.best_AR_stage[3].p == 3
    @test par_1_fixed_p.best_AR_stage[4].p == 3
    @test par_1_fixed_p.best_AR_stage[5].p == 3
    @test par_1_fixed_p.best_AR_stage[6].p == 3
    @test par_1_fixed_p.best_AR_stage[7].p == 3
    @test par_1_fixed_p.best_AR_stage[8].p == 3
    @test par_1_fixed_p.best_AR_stage[9].p == 3
    @test par_1_fixed_p.best_AR_stage[10].p == 3
    @test par_1_fixed_p.best_AR_stage[11].p == 3
    @test par_1_fixed_p.best_AR_stage[12].p == 3

     # TODO test something about the simulations
    scen_b, scen_f = PAR.simulate_par_f_b([par_1], 10, 100, 20)
    scen_b, scen_f = PAR.simulate_par_f_b([par_1; par_2], 10, 100, 20)
    scen_b, scen_f = PAR.simulate_par_f_b([par_1; par_2], 10, 100, 20; global_lower_bound = true)
end