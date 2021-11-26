#! format: off
# This file is gonna be a hot mess anyway, no reason to format it
using Test, CSV

include("utils.jl")

@testset "fund_args" begin
    @test @testbench fund_args(t) (t::Range{Float64,0,1e7})
end
@testset "norm_ang" begin
    @test @testbench norm_ang(θ) (θ::Range{Float64,-1e5,1e5})
end
@testset "ee_ct" begin
    @test @testbench ee_ct(jd_high, jd_low; accuracy = :full) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
    @test @testbench ee_ct(jd_high, jd_low; accuracy = :reduced) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
end

@testset "nu2000k" begin
    @test @testbench nu2000k(jd_high, jd_low) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
end
@testset "iau2000a" begin
    @test @testbench iau2000a(jd_high, jd_low) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
end

@testset "nutation_angles" begin
    @test @testbench nutation_angles(t; accuracy = :full) (t::Range{Float64,-1e3,4e3})
    @test @testbench nutation_angles(t; accuracy = :reduced) (t::Range{Float64,-1e3,4e3})
end
@testset "mean_obliq" begin
    @test @testbench mean_obliq(jd_tdb) (jd_tdb::Range{Float64,0,1e7})
end
@testset "frame_tie" begin
    @test @testbench frame_tie(pos, :dynamic2icrs) (pos::Position)
    @test @testbench frame_tie(pos, :icrs2dynamic) (pos::Position)
end
@testset "nutation" begin
    @test @testbench nutation(jd_tdb, pos; accuracy = :full, direction = :mean2true) (jd_tdb::Range{Float64,-1e7,1e7}) (pos::Position)
    @test @testbench nutation(jd_tdb, pos; accuracy = :reduced, direction = :mean2true) (jd_tdb::Range{Float64,-1e7,1e7}) (pos::Position)
    @test @testbench nutation(jd_tdb, pos; accuracy = :full, direction = :true2mean) (jd_tdb::Range{Float64,-1e7,1e7}) (pos::Position)
    @test @testbench nutation(jd_tdb, pos; accuracy = :reduced, direction = :true2mean) (jd_tdb::Range{Float64,-1e7,1e7}) (pos::Position)
end
@testset "precession" begin
    @test @testbench precession(jd_tdb, pos, NOVAS.T0) (jd_tdb::Range{Float64,-1e7,1e7}) (pos::Position)
end
@testset "e_tilt" begin
    @test @testbench e_tilt(jd_tdb; accuracy = :full) (jd_tdb::Range{Float64,-1e7,1e7})
    @test @testbench e_tilt(jd_tdb; accuracy = :reduced) (jd_tdb::Range{Float64,-1e7,1e7})
end
@testset "ira_equinox" begin
    @test @testbench ira_equinox(jd_tdb; accuracy = :full, equinox = :mean) (jd_tdb::Range{Float64,-1e7,1e7})
    @test @testbench ira_equinox(jd_tdb; accuracy = :full, equinox = :true) (jd_tdb::Range{Float64,-1e7,1e7})
    @test @testbench ira_equinox(jd_tdb; accuracy = :reduced, equinox = :mean) (jd_tdb::Range{Float64,-1e7,1e7})
    @test @testbench ira_equinox(jd_tdb; accuracy = :reduced, equinox = :true) (jd_tdb::Range{Float64,-1e7,1e7})
end
@testset "cio_location" begin
    @test @testbench cio_location(jd_tdb; accuracy = :full) (jd_tdb::Range{Float64,-1e7,1e7})
end
@testset "cio_basis" begin
    @test @testbench cio_basis(jd_tdb, ra_cio; accuracy = :full) (jd_tdb::Range{Float64,-1e7,1e7}) (ra_cio::Range{Float64,0,24})
    @test @testbench cio_basis(jd_tdb, ra_cio; accuracy = :reduced) (jd_tdb::Range{Float64,-1e7,1e7}) (ra_cio::Range{Float64,0,24})
end
@testset "era" begin
    @test @testbench era(jd_high, jd_low) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
end
@testset "tdb2tt" begin
    @test @testbench tdb2tt(jdb_jd) (jdb_jd::Range{Float64,-1e7,1e7})
end
@testset "sidereal_time" begin
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :mean, method = :CIO, accuracy = :full) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :mean, method = :CIO, accuracy = :reduced) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :mean, method = :equinox, accuracy = :full) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :mean, method = :equinox, accuracy = :reduced) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :apparent, method = :CIO, accuracy = :full) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :apparent, method = :CIO, accuracy = :reduced) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :apparent, method = :equinox, accuracy = :full) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
    @test @testbench sidereal_time(jd_high, jd_low, delta_t; gst_type = :apparent, method = :equinox, accuracy = :reduced) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70})
end
@testset "wobble" begin
    @test @testbench wobble(tjd, xp, yp, pos; direction = :itrs2terr) (tjd::Range{Float64,-1e7,1e7}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
    @test @testbench wobble(tjd, xp, yp, pos; direction = :terr2itrs) (tjd::Range{Float64,-1e7,1e7}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
end
@testset "spin" begin
    @test @testbench spin(θ, pos) (θ::Range{Float64,0,360}) (pos::Position)
end
@testset "ter2cel" begin
    @test @testbench ter2cel(jd_high, jd_low, delta_t, pos; method = :CIO, accuracy = :full, option = :GCRS, xp = xp, yp = yp) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
    @test @testbench ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :full, option = :GCRS, xp = xp, yp = yp) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
    @test @testbench ter2cel(jd_high, jd_low, delta_t, pos; method = :CIO, accuracy = :reduced, option = :GCRS, xp = xp, yp = yp) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
    @test @testbench ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :reduced, option = :GCRS, xp = xp, yp = yp) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
    # Option only matters for method = :equinox
    @test @testbench ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :full, option = :equinox, xp = xp, yp = yp) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
    @test @testbench ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :reduced, option = :equinox, xp = xp, yp = yp) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1}) (delta_t::Range{Float64,-10,70}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1}) (pos::Position)
end
@testset "refract" begin
    @test @testbench refract(location, zd; ref_option = :standard) (location::Location) (zd::Range{Float64,0,90})
    @test @testbench refract(location, zd; ref_option = :location) (location::Location) (zd::Range{Float64,0,90})
end
@testset "equ2hor" begin
    @test @testbench equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy = :full, xp = xp, yp = yp, ref_option = :none) (jd_ut1::Range{Float64,0,1e7}) (delta_t::Range{Float64,-10,70}) (ra::Range{Float64,0,24}) (dec::Range{Float64,0,360}) (location::Location) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
    @test @testbench equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy = :full, xp = xp, yp = yp, ref_option = :standard) (jd_ut1::Range{Float64,0,1e7}) (delta_t::Range{Float64,-10,70}) (ra::Range{Float64,0,24}) (dec::Range{Float64,0,360}) (location::Location) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
    @test @testbench equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy = :full, xp = xp, yp = yp, ref_option = :location) (jd_ut1::Range{Float64,0,1e7}) (delta_t::Range{Float64,-10,70}) (ra::Range{Float64,0,24}) (dec::Range{Float64,0,360}) (location::Location) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
    @test @testbench equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy = :reduced, xp = xp, yp = yp, ref_option = :none) (jd_ut1::Range{Float64,0,1e7}) (delta_t::Range{Float64,-10,70}) (ra::Range{Float64,0,24}) (dec::Range{Float64,0,360}) (location::Location) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
    @test @testbench equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy = :reduced, xp = xp, yp = yp, ref_option = :standard) (jd_ut1::Range{Float64,0,1e7}) (delta_t::Range{Float64,-10,70}) (ra::Range{Float64,0,24}) (dec::Range{Float64,0,360}) (location::Location) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
    @test @testbench equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy = :reduced, xp = xp, yp = yp, ref_option = :location) (jd_ut1::Range{Float64,0,1e7}) (delta_t::Range{Float64,-10,70}) (ra::Range{Float64,0,24}) (dec::Range{Float64,0,360}) (location::Location) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
end

# Serialize Results
CSV.write("../benchmarks.csv", results)
