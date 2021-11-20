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
    @testset "mean_obliq" begin
        @test mean_obliq(jd_high) ≈ NOVAS.mean_obliq(jd_high)
    end
    # Generate random position vector
    pos = rand(3)
    @testset "frame_tie" begin
        @test frame_tie(pos,-1) ≈ NOVAS.frame_tie(pos,:dynamic2icrs)
        @test frame_tie(pos,0) ≈ NOVAS.frame_tie(pos,:icrs2dynamic)
    end
    @testset "nutation" begin
        @test nutation(jd_high,0,0,pos) ≈ NOVAS.nutation(jd_high,pos;accuracy=:full,direction=:mean2true)
        @test nutation(jd_high,0,1,pos) ≈ NOVAS.nutation(jd_high,pos;accuracy=:reduced,direction=:mean2true)
        @test nutation(jd_high,1,0,pos) ≈ NOVAS.nutation(jd_high,pos;accuracy=:full,direction=:true2mean)
        @test nutation(jd_high,1,1,pos) ≈ NOVAS.nutation(jd_high,pos;accuracy=:reduced,direction=:true2mean)
    end
    @testset "precession" begin
        @test precession(jd_high,pos,NOVAS.T0) ≈ NOVAS.precession(jd_high,pos,NOVAS.T0)
        @test precession(NOVAS.T0,pos,jd_high) ≈ NOVAS.precession(NOVAS.T0,pos,jd_high)
    end
    @testset "e_tilt" begin
        @test e_tilt(jd_high,0) ≈ NOVAS.e_tilt(jd_high)
        @test e_tilt(jd_high,1) ≈ NOVAS.e_tilt(jd_high;accuracy=:reduced)
    end
    @testset "ira_equinox" begin
        @test ira_equinox(jd_high,0,0) ≈ NOVAS.ira_equinox(jd_high;accuracy=:full,equinox=:mean)
        @test ira_equinox(jd_high,0,1) ≈ NOVAS.ira_equinox(jd_high;accuracy=:reduced,equinox=:mean)
        @test ira_equinox(jd_high,1,0) ≈ NOVAS.ira_equinox(jd_high;accuracy=:full,equinox=:true)
        @test ira_equinox(jd_high,1,1) ≈ NOVAS.ira_equinox(jd_high;accuracy=:reduced,equinox=:true)
    end
    @testset "cio_location" begin
        @test cio_location(jd_high,0) ≈ NOVAS.cio_location(jd_high;accuracy=:full)
        @test cio_location(jd_high,1) ≈ NOVAS.cio_location(jd_high;accuracy=:reduced)
    end
    # Random RA in hours
    ra_cio = rand() * 24
    @testset "cio_basis" begin
        @test cio_basis(jd_high,ra_cio,2,0) ≈ NOVAS.cio_basis(jd_high,ra_cio;accuracy=:full)
        @test cio_basis(jd_high,ra_cio,2,1) ≈ NOVAS.cio_basis(jd_high,ra_cio;accuracy=:reduced)
    end
    @testset "era" begin
        @test era(jd_high,jd_low) ≈ NOVAS.era(jd_high,jd_low)
    end
    @testset "tdb2tt" begin
        @test tdb2tt(jd_high) ≈ NOVAS.tdb2tt(jd_high)
    end
    @testset "sidereal_time" begin
        @test sidereal_time(jd_high,jd_low,0.0,0,0,0) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:mean,method=:CIO,accuracy=:full)
        @test sidereal_time(jd_high,jd_low,0.0,0,0,1) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:mean,method=:CIO,accuracy=:reduced)
        @test sidereal_time(jd_high,jd_low,0.0,0,1,0) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:mean,method=:equinox,accuracy=:full)
        @test sidereal_time(jd_high,jd_low,0.0,0,1,1) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:mean,method=:equinox,accuracy=:reduced)
        @test sidereal_time(jd_high,jd_low,0.0,1,0,0) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:apparent,method=:CIO,accuracy=:full)
        @test sidereal_time(jd_high,jd_low,0.0,1,0,1) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:apparent,method=:CIO,accuracy=:reduced)
        @test sidereal_time(jd_high,jd_low,0.0,1,1,0) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:apparent,method=:equinox,accuracy=:full)
        @test sidereal_time(jd_high,jd_low,0.0,1,1,1) ≈ NOVAS.sidereal_time(jd_high,jd_low,0.0;gst_type=:apparent,method=:equinox,accuracy=:reduced)
    end
end