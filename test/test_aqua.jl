module TestAqua

using Aqua
using PeriodicAutoregressive
using Test

function runtests()
    @testset "Aqua" begin
        @testset "Ambiguities" begin
            Aqua.test_ambiguities(PeriodicAutoregressive, recursive = false)
        end
        Aqua.test_all(PeriodicAutoregressive, ambiguities = false)
    end
end

TestAqua.runtests()

end