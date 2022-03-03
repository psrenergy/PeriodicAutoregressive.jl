@testset "PAR(p)-A" begin
    n_stages = 12
    p_lim = 6
    par_p_a_1 = PARpA(funil_grande, n_stages, p_lim; information_criteria = "aic");
    par_p_a_2 = PARpA(camargos, n_stages, p_lim; information_criteria = "aic");
    fit_par!(par_p_a_1)
    fit_par!(par_p_a_2)

    # TODO get all the estima files in this repository
    # Some random tests based on an estima.csv file
    @test par_p_a_1.μ_stage[2] ≈ 286.75 atol = 1e-2
    @test par_p_a_1.σ_stage[2] ≈ 124.45 atol = 1e-2
    @test par_p_a_1.best_AR_A_stage[4].p == 2
    @test par_p_a_2.y_anual[13] ≈ 233.1666666
    @test par_p_a_2.y_anual[14] == 255.75
    @test par_p_a_2.y_anual[end] == 91.25
    
    # TODO test something about the simulations
    scen = simulate_par(par_p_a_1, 10, 100)
    scen = simulate_par([par_p_a_1; par_p_a_2], 100, 1000)

    p_lim = 3
    par_1_fixed_p = PARpA(funil_grande, n_stages, p_lim; information_criteria = "fixed_at_p_lim");
    fit_par!(par_1_fixed_p);
    @test par_1_fixed_p.best_AR_A_stage[1].p == 3
    @test par_1_fixed_p.best_AR_A_stage[2].p == 3
    @test par_1_fixed_p.best_AR_A_stage[3].p == 3
    @test par_1_fixed_p.best_AR_A_stage[4].p == 3
    @test par_1_fixed_p.best_AR_A_stage[5].p == 3
    @test par_1_fixed_p.best_AR_A_stage[6].p == 3
    @test par_1_fixed_p.best_AR_A_stage[7].p == 3
    @test par_1_fixed_p.best_AR_A_stage[8].p == 3
    @test par_1_fixed_p.best_AR_A_stage[9].p == 3
    @test par_1_fixed_p.best_AR_A_stage[10].p == 3
    @test par_1_fixed_p.best_AR_A_stage[11].p == 3
    @test par_1_fixed_p.best_AR_A_stage[12].p == 3

    p_lim = 20
    @test_throws AssertionError PARpA(funil_grande, n_stages, p_lim; information_criteria = "fixed_at_p_lim");
end