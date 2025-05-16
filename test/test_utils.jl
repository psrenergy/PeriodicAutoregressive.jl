module TestUtils

using PeriodicAutoregressive
using Test

function test_assert_series_without_missing()
    @test_throws AssertionError PeriodicAutoregressive.assert_series_without_missing([1, 2, NaN])
    @test PeriodicAutoregressive.assert_series_without_missing([1.0, 2, 3])
end

function test_vov()
    vov = [[1, 2, 3, 4], [1, 2]]
    m = PeriodicAutoregressive.concatenate_from_the_bottom_elements(vov)
    @test m == [3 1; 4 2]
end

function test_concatenate_from_the_bottom_elements()
    vov = [rand(15), rand(10), rand(100)]
    m = PeriodicAutoregressive.concatenate_from_the_bottom_elements(vov)
    @test size(m) == (10, 3)
end

function runtests()
    Base.GC.gc()
    Base.GC.gc()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

TestUtils.runtests()

end