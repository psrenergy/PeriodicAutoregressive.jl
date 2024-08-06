@testset "PAR(p)" begin
    n_stages = 12
    p_lim = 6
    par_1 = PARp(funil_grande, n_stages, p_lim; information_criteria = "aic");
    par_2 = PARp(batalha, n_stages, p_lim; information_criteria = "aicc");
    par_3 = PARp(zeros(length(funil_grande)), n_stages, p_lim);
    par_4 = PARp(inflow_with_stages_0_variance, n_stages, p_lim; information_criteria = "aic");
    fit_par!(par_1);
    fit_par!(par_2);
    fit_par!(par_3);
    fit_par!(par_4);

    @test par_1.μ_stage[2] ≈ 286.75 atol = 1e-2
    @test par_1.σ_stage[2] ≈ 124.45 atol = 1e-2
    @test par_1.best_AR_stage[4].p == 4

    @test par_2.best_AR_stage[2].ϕ[1] .≈ 0.66114 atol = 1e-2
    @test par_2.best_AR_stage[3].ϕ[1] .≈ 0.45045 atol = 1e-2

    @test par_3.μ_stage[1] == 0.0
    @test par_3.σ_stage[2] == 0.0

    scen = simulate_par(par_1, 10, 100)
    scen = simulate_par([par_1, par_2], 100, 1000)
    scen = simulate_par([par_1, par_2, par_3], 100, 1000)
    @test all(scen[:, 3, :] .== 0.0)

    scen = simulate_par([par_1, par_3, par_2], 100, 1000)
    @test all(scen[:, 2, :] .== 0.0)

    scen = simulate_par([par_3], 10, 12)
    scen_b, lags = PAR.simulate_par_f_b([par_3], 10, 10, 12)

    @test sum(scen) == 0
    @test sum(scen_b) == 0

    fit_par!(par_4);
    scen = simulate_par([par_4], 10, 12)
    scen_b, lags = PAR.simulate_par_f_b([par_4], 10, 10, 12)

    @test isempty(findall(isnan, scen))
    @test isempty(findall(isnan, scen_b))

    scen = simulate_par([par_1, par_2, par_3, par_4], 100, 1000)
    @test all(scen[7, 4, :] .== 1100.0)
    @test all(scen[:, 3, :] .== 0.0)
    scen = simulate_par([par_1, par_4, par_2, par_3], 100, 1000)
    @test all(scen[7, 2, :] .== 1100.0)
    @test all(scen[:, 4, :] .== 0.0)
    scen = simulate_par([par_4, par_3, par_2, par_1], 100, 1000)
    @test all(scen[7, 1, :] .== 1100.0)
    @test all(scen[:, 2, :] .== 0.0)

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

    scen_b, scen_f = PAR.simulate_par_f_b([par_1], 10, 100, 20)
    scen_b, scen_f = PAR.simulate_par_f_b([par_1; par_2], 10, 100, 20)
    scen_b, scen_f = PAR.simulate_par_f_b([par_1; par_2], 10, 100, 20; has_global_lower_bound = true)
end