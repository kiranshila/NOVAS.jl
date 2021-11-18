using NOVAS
using Test

# Add custom isapprox to test tuples
import Base.isapprox
Base.isapprox(x::Tuple, y::Tuple; kws...) = isapprox(collect(x), collect(y); kws...)

# Include un-exported c wrappers
include("wrapper.jl")

@testset "utilities" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    # Figure the TDB Julian centuries
    t = rand()
    @test fund_args(t) ≈ NOVAS.fund_args(t)
    # Generate random big angle
    θ = rand() * 100
    @test norm_ang(θ) ≈ NOVAS.norm_ang(θ)
    @test norm_ang(-θ) ≈ NOVAS.norm_ang(-θ)
    # Test ee_ct
    @test ee_ct(jd_high,jd_low,0) ≈ NOVAS.ee_ct(jd_high,jd_low)
    @test ee_ct(jd_high,jd_low,1) ≈ NOVAS.ee_ct(jd_high,jd_low;accuracy=:reduced)
end

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
