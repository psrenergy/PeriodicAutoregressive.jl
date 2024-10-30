function test_aqua()
    @testset "Ambiguities" begin
        Aqua.test_ambiguities(PeriodicAutoregressive, recursive = false)
    end
    Aqua.test_all(PeriodicAutoregressive, ambiguities = false)

    return nothing
end
