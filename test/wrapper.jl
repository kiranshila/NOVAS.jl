using NOVAS
import NOVAS.LibNOVAS

function c_fund_args(t)
    a = zeros(Cdouble, 5)
    LibNOVAS.fund_args(t, a)
    return a
end

function norm_ang(angle)
    LibNOVAS.norm_ang(angle)
end

function ee_ct(jd_high, jd_low, accuracy)
    LibNOVAS.ee_ct(jd_high, jd_low, accuracy)
end

function iau2000a(jd_high, jd_low)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.iau2000a(jd_high, jd_low, dpsi, deps)
    return dpsi[], deps[]
end

function nu2000k(jd_high, jd_low)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.nu2000k(jd_high, jd_low, dpsi, deps)
    return dpsi[], deps[]
end

# Accuracy flags, 0 - Full, 1 - Reduced

function nutation_angles(t, accuracy)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.nutation_angles(t, accuracy, dpsi, deps)
    return dpsi[], deps[]
end

function e_tilt(jd_tdb, accuracy)
    mobl = Ref{Cdouble}(0.0)
    tobl = Ref{Cdouble}(0.0)
    ee = Ref{Cdouble}(0.0)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.e_tilt(jd_tdb, accuracy, mobl, tobl, ee, dpsi, deps)
    return mobl[], tobl[], ee[], dpsi[], deps[]
end

function mean_obliq(jd_tdb)
    LibNOVAS.mean_obliq(jd_tdb)
end

function frame_tie(pos, direction)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.frame_tie(pos, direction, pos2)
    return pos2
end

function nutation(jd_tdb, direction, accuracy, pos)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.nutation(jd_tdb, direction, accuracy, pos, pos2)
    return pos2
end

function precession(jd_tdb1, pos, jd_tdb2)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.precession(jd_tdb1, pos, jd_tdb2, pos2)
    return pos2
end

function ira_equinox(jd_tdb, equinox, accuracy)
    LibNOVAS.ira_equinox(jd_tdb, equinox, accuracy)
end

function cio_location(jd_tdb, accuracy)
    ra_cio = Ref{Cdouble}(0.0)
    ref_sys = Ref{Cshort}(0)
    LibNOVAS.cio_location(jd_tdb, accuracy, ra_cio, ref_sys)
    return ra_cio[]
end

function cio_basis(jd_tdb, ra_cio, ref_sys, accuracy)
    x = zeros(Cdouble, 3)
    y = zeros(Cdouble, 3)
    z = zeros(Cdouble, 3)
    LibNOVAS.cio_basis(jd_tdb, ra_cio, ref_sys, accuracy, x, y, z)
    return x, y, z
end

function era(jd_high, jd_low)
    LibNOVAS.era(jd_high, jd_low)
end

function tdb2tt(jdb_jd)
    tt_jd = Ref{Cdouble}(0.0)
    secdiff = Ref{Cdouble}(0.0)
    LibNOVAS.tdb2tt(jdb_jd, tt_jd, secdiff)
    return tt_jd[], secdiff[]
end

function sidereal_time(jd_high,jd_low,delta_t,gst_type,method,accuracy)
    gst = Ref{Cdouble}(0.0)
    LibNOVAS.sidereal_time(jd_high, jd_low, delta_t, gst_type, method, accuracy, gst)
    return gst[]
end

function wobble(tjd, direction, xp, yp, pos)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.wobble(tjd, direction, xp, yp, pos, pos2)
    return pos2
end

function spin(angle, pos)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.spin(angle, pos, pos2)
    return pos2

end

function ter2cel(jd_ut_high,jd_ut_low,delta_t,method,accuracy,option,xp,yp,vec)
    vec2 = zeros(Cdouble, 3)
    LibNOVAS.ter2cel(jd_ut_high, jd_ut_low, delta_t, method, accuracy, option, xp, yp, vec, vec2)
    return vec2
end

function refract(location, ref_option::Int, zd_obs::Real)
    LibNOVAS.refract(location, ref_option, zd_obs)
end

function equ2hor(jd_ut1,delta_t,accuracy,xp,yp,location,ra,dec,ref_option)
    zd = Ref{Cdouble}(0.0)
    az = Ref{Cdouble}(0.0)
    rar = Ref{Cdouble}(0.0)
    decr = Ref{Cdouble}(0.0)
    LibNOVAS.equ2hor(jd_ut1, delta_t, accuracy, xp, yp, location, ra, dec, ref_option, zd, az, rar, decr)
    return zd[], az[], rar[], decr[]
end

function test_equ2hor(accuracy,ref_option,xp,yp)
    if accuracy == 0
        acc = :full
    elseif accuracy == 1
        acc = :reduced
    end

    if ref_option == 0
        ref = :none
    elseif ref_option == 1
        ref = :standard
    elseif ref_option == 2
        ref = :location
    end

    # Random location
    lat = rand(Cdouble) * 180 - 90
    lon = rand(Cdouble) * 180 - 90
    alt = rand(Cdouble) * 8000
    temp = rand(Cdouble) * 50
    pressure = rand(Cdouble) * 1000
    location  = NOVAS.OnSurface(lat, lon, alt, temp, pressure)
    clocation = Ref{LibNOVAS.on_surface}(LibNOVAS.on_surface(lat, lon, alt, temp, pressure))

    # Rand jd_high
    jd_high = rand()*1e6
    delta_t = rand()

    # Random sky position
    ra = rand() * 24
    dec = rand() * 180 - 90

    c = equ2hor(jd_high, delta_t, accuracy, xp, yp, clocation, ra, dec, ref_option)
    julia = NOVAS.equ2hor(jd_high, delta_t, ra, dec, location; accuracy = acc, xp = xp, yp = yp, ref_option = ref)
    return c ≈ julia
end