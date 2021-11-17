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
    @test nu2000k(jd_high,jd_low) ≈ NOVAS.nu2000k(jd_high,jd_low)
    @test iau2000a(jd_high,jd_low) ≈ NOVAS.iau2000a(jd_high,jd_low)
end

@testset "novas" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    # Figure the TDB Julian centuries
    t = ((jd_high - NOVAS.T0) + jd_low) / 36525.0
    @test nutation_angles(t,0) ≈ NOVAS.nutation_angles(t)
    @test nutation_angles(t,1) ≈ NOVAS.nutation_angles(t;accuracy=:reduced)
end
