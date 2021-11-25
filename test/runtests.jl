using Test, CSV

include("utils.jl")

@testset verbose = true "utilities" begin
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
end

@testset verbose = true "nutation" begin
    # Generate a julian date
    @testset "nu2000k" begin
        @test @testbench nu2000k(jd_high, jd_low) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
    end
    @testset "iau2000a" begin
        @test @testbench iau2000a(jd_high, jd_low) (jd_high::Range{Float64,0,1e7}) (jd_low::Range{Float64,0,1})
    end
end

@testset verbose = true "novas" begin
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
        @test @testbench wobble(tjd,xp, yp, pos; direction = :itrs2terr) (tjd::Range{Float64,-1e7,1e7}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
        @test @testbench wobble(tjd,xp, yp, pos; direction = :terr2itrs) (tjd::Range{Float64,-1e7,1e7}) (xp::Range{Float64,-1,1}) (yp::Range{Float64,-1,1})
    end
    #=
    # Random angle
    angle = rand()
    @testset "spin" begin
        @test spin(angle, pos) ≈ NOVAS.spin(angle, pos)
    end
    @testset "ter2cel" begin
        @test ter2cel(jd_high, jd_low, delta_t, 0, 0, 0, xp, yp, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :CIO, accuracy = :full, option = :GCRS, xp = xp, yp = yp)
        @test ter2cel(jd_high, jd_low, delta_t, 1, 0, 0, xp, yp, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :full, option = :GCRS, xp = xp, yp = yp)
        @test ter2cel(jd_high, jd_low, delta_t, 1, 0, 1, xp, yp, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :full, option = :equinox, xp = xp, yp = yp)
        @test ter2cel(jd_high, jd_low, delta_t, 0, 1, 0, xp, yp, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :CIO, accuracy = :reduced, option = :GCRS, xp = xp, yp = yp)
        @test ter2cel(jd_high, jd_low, delta_t, 1, 1, 0, xp, yp, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :reduced, option = :GCRS, xp = xp, yp = yp)
        @test ter2cel(jd_high, jd_low, delta_t, 1, 1, 1, xp, yp, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :reduced, option = :equinox, xp = xp, yp = yp)
        # Now with zero xp and yp
        @test ter2cel(jd_high, jd_low, delta_t, 1, 1, 0, 0.0, 0.0, pos) ≈ NOVAS.ter2cel(jd_high, jd_low, delta_t, pos; method = :equinox, accuracy = :reduced, option = :GCRS)
    end
    # Generate a random location
    lat = rand(Cdouble) * 180 - 90
    lon = rand(Cdouble) * 180 - 90
    alt = rand(Cdouble) * 8000
    temp = rand(Cdouble) * 50
    pressure = rand(Cdouble) * 1000
    location = NOVAS.OnSurface(lat, lon, alt, temp, pressure)
    clocation = Ref{LibNOVAS.on_surface}(LibNOVAS.on_surface(lat, lon, alt, temp, pressure))
    zd = rand() * 90
    @testset "refract" begin
        @test refract(clocation, 1, zd) ≈ NOVAS.refract(location, zd; ref_option = :standard)
        @test refract(clocation, 2, zd) ≈ NOVAS.refract(location, zd; ref_option = :location)
    end
    # Random coordinate
    ra = rand() * 24
    dec = rand() * 180 - 90
    @testset "equ2hor" begin
        @test equ2hor(jd_high, delta_t, 0, xp, yp, clocation, ra, dec, 0) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :full, xp = xp, yp = yp, ref_option = :none)
        @test equ2hor(jd_high, delta_t, 0, xp, yp, clocation, ra, dec, 1) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :full, xp = xp, yp = yp, ref_option = :standard)
        @test equ2hor(jd_high, delta_t, 0, xp, yp, clocation, ra, dec, 2) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :full, xp = xp, yp = yp, ref_option = :location)
        @test equ2hor(jd_high, delta_t, 1, xp, yp, clocation, ra, dec, 0) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :reduced, xp = xp, yp = yp, ref_option = :none)
        @test equ2hor(jd_high, delta_t, 1, xp, yp, clocation, ra, dec, 1) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :reduced, xp = xp, yp = yp, ref_option = :standard)
        @test equ2hor(jd_high, delta_t, 1, xp, yp, clocation, ra, dec, 2) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :reduced, xp = xp, yp = yp, ref_option = :location)
        # And with zero polar motion
        @test equ2hor(jd_high, delta_t, 0, 0.0, 0.0, clocation, ra, dec, 2) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :full, ref_option = :location)
        @test equ2hor(jd_high, delta_t, 1, 0.0, 0.0, clocation, ra, dec, 2) ≈ NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = :reduced, ref_option = :location)
    end
     =#
end


# Serialize Results
CSV.write("benchmarks.csv", results)