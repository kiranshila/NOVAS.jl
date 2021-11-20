using NOVAS
using Test

# Add custom isapprox to test tuples
import Base.isapprox
Base.isapprox(x::Tuple, y::Tuple; kws...) = isapprox(collect(x), collect(y); kws...)

# Include un-exported c wrappers
include("wrapper.jl")

@testset verbose = true "utilities" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    # Figure the TDB Julian centuries
    t = rand()
    # Generate random big angle
    θ = rand() * 100
    @testset "fund_args" begin
        @test fund_args(t) ≈ NOVAS.fund_args(t)
    end
    @testset "norm_ang" begin
        @test norm_ang(θ) ≈ NOVAS.norm_ang(θ)
        @test norm_ang(-θ) ≈ NOVAS.norm_ang(-θ)
    end
    @testset "ee_ct" begin
        @test ee_ct(jd_high,jd_low,0) ≈ NOVAS.ee_ct(jd_high,jd_low)
        @test ee_ct(jd_high,jd_low,1) ≈ NOVAS.ee_ct(jd_high,jd_low;accuracy=:reduced)
    end
end

@testset verbose = true "nutation" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    @testset "nu2000k" begin
        @test nu2000k(jd_high,jd_low) ≈ NOVAS.nu2000k(jd_high,jd_low)
    end
    @testset "iau2000a" begin
        @test iau2000a(jd_high,jd_low) ≈ NOVAS.iau2000a(jd_high,jd_low)
    end
end

@testset verbose = true "novas" begin
    # Generate a julian date
    jd_high = rand()*1e6 |> round
    jd_low = rand()
    # Figure the TDB Julian centuries
    t = ((jd_high - NOVAS.T0) + jd_low) / 36525.0
    @testset "nutation_angles" begin
        @test nutation_angles(t,0) ≈ NOVAS.nutation_angles(t)
        @test nutation_angles(t,1) ≈ NOVAS.nutation_angles(t;accuracy=:reduced)
    end
    @testset "e_tilt" begin
        @test e_tilt(jd_high,0) ≈ NOVAS.e_tilt(jd_high)
        @test e_tilt(jd_high,1) ≈ NOVAS.e_tilt(jd_high;accuracy=:reduced)
    end
    @testset "mean_obliq" begin
        @test mean_obliq(jd_high) ≈ NOVAS.mean_obliq(jd_high)
    end
    # Generate random position vector
    pos = rand(3)
    @testset "frame_tie" begin
        @test frame_tie(pos,-1) ≈ NOVAS.frame_tie(pos,:dynamic2icrs)
        @test frame_tie(pos,0) ≈ NOVAS.frame_tie(pos,:icrs2dynamic)
    end
end
