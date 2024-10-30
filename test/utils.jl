function test_utils()
    @test_throws AssertionError PeriodicAutoregressive.assert_series_without_missing([1,2,NaN])
    @test PeriodicAutoregressive.assert_series_without_missing([1.0,2,3])

    vov = [[1,2,3,4] , [1,2]]
    m = PeriodicAutoregressive.concatenate_from_the_bottom_elements(vov)
    @test m == [3 1; 4 2]

    vov = [rand(15), rand(10), rand(100)]
    m = PeriodicAutoregressive.concatenate_from_the_bottom_elements(vov)
    @test size(m) == (10, 3)

    return nothing
end