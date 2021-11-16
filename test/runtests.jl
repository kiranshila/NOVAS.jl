using NOVAS
using Test

# Add custom isapprox to test tuples
import Base.isapprox
Base.isapprox(x::Tuple, y::Tuple; kws...) = isapprox(collect(x), collect(y); kws...)

# Include un-exported c wrappers
include("wrapper.jl")

@testset "nutation" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    # Test nutation conversions
    @test iau2000a(jd_high,jd_low) ≈ NOVAS.iau2000a(jd_high,jd_low) atol=eps(Float64)
end

@testset "novas" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    # Figure the TDB Julian centuries
    t = ((jd_high - NOVAS.T0) + jd_low) / 36525.0
end
