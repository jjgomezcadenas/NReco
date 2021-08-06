using NReco
using Test

@testset "NReco.jl" begin
    x = collect(1:100)
    y = collect(6:9)
    xr = NReco.in_range(x, 5, 10) # interval is [ )
    println(xr)
    println(y)
    @test all(y .== xr)

end
