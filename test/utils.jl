@testset "utils.jl" begin
    @test_throws AssertionError PAR.assert_series_without_missing([1,2,NaN])
    @test PAR.assert_series_without_missing([1.0,2,3])

    vov = [[1,2,3,4] , [1,2]]
    m = PAR.concatenate_from_the_bottom_elements(vov)
    @test m == [3 1; 4 2]

    vov = [rand(15), rand(10), rand(100)]
    m = PAR.concatenate_from_the_bottom_elements(vov)
    @test size(m) == (10, 3)
end