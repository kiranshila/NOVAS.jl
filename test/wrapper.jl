using NOVAS
import NOVAS.LibNOVAS

function c_fund_args(t::Real)
    a = zeros(Cdouble, 5)
    LibNOVAS.fund_args(t, a)
    return a
end

function norm_ang(angle::Real)
    LibNOVAS.norm_ang(angle)
end

function ee_ct(jd_high::Real, jd_low::Real, accuracy::Int)
    LibNOVAS.ee_ct(jd_high, jd_low, accuracy)
end

function iau2000a(jd_high::Real, jd_low::Real)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.iau2000a(jd_high, jd_low, dpsi, deps)
    return dpsi[], deps[]
end

function nu2000k(jd_high::Real, jd_low::Real)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.nu2000k(jd_high, jd_low, dpsi, deps)
    return dpsi[], deps[]
end

# Accuracy flags, 0 - Full, 1 - Reduced

function nutation_angles(t::Real, accuracy::Int)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.nutation_angles(t, accuracy, dpsi, deps)
    return dpsi[], deps[]
end

function e_tilt(jd_tdb::Real, accuracy::Int)
    mobl = Ref{Cdouble}(0.0)
    tobl = Ref{Cdouble}(0.0)
    ee = Ref{Cdouble}(0.0)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.e_tilt(jd_tdb, accuracy, mobl, tobl, ee, dpsi, deps)
    return mobl[], tobl[], ee[], dpsi[], deps[]
end

function mean_obliq(jd_tdb::Real)
    LibNOVAS.mean_obliq(jd_tdb)
end

function frame_tie(pos::AbstractVector, direction::Int)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.frame_tie(pos, direction, pos2)
    return pos2
end

function nutation(jd_tdb::Real, direction::Int, accuracy::Int, pos::AbstractVector)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.nutation(jd_tdb, direction, accuracy, pos, pos2)
    return pos2
end

function precession(jd_tdb1::Real, pos::AbstractVector, jd_tdb2::Real)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.precession(jd_tdb1, pos, jd_tdb2, pos2)
    return pos2
end

function ira_equinox(jd_tdb::Real, equinox::Int, accuracy::Int)
    LibNOVAS.ira_equinox(jd_tdb, equinox, accuracy)
end

function cio_location(jd_tdb::Real, accuracy::Int)
    ra_cio = Ref{Cdouble}(0.0)
    ref_sys = Ref{Cdouble}(0.0)
    LibNOVAS.cio_location(jd_tdb, accuracy, ra_cio, ref_sys)
    return ra_cio[]
end

function cio_basis(jd_tdb::Real, ra_cio::Real, ref_sys::Int, accuracy::Int)
    x = zeros(Cdouble, 3)
    y = zeros(Cdouble, 3)
    z = zeros(Cdouble, 3)
    LibNOVAS.cio_basis(jd_tdb, ra_cio, ref_sys, accuracy, x, y, z)
    return x, y, z
end

function era(jd_high::Real, jd_low::Real)
    LibNOVAS.era(jd_high, jd_low)
end

function tdb2tt(jdb_jd::Real)
    tt_jd = Ref{Cdouble}(0.0)
    secdiff = Ref{Cdouble}(0.0)
    LibNOVAS.tdb2tt(jdb_jd, tt_jd, secdiff)
    return tt_jd[], secdiff[]
end

function sidereal_time(jd_high::Real,
    jd_low::Real,
    delta_t::Real,
    gst_type::Int,
    method::Int,
    accuracy::Int)
    gst = Ref{Cdouble}(0.0)
    LibNOVAS.sidereal_time(jd_high, jd_low, delta_t, gst_type, method, accuracy, gst)
    return gst[]
end

function wobble(tjd::Real, direction::Int, xp::Real, yp::Real, pos::Vector)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.wobble(tjd, direction, xp, yp, pos, pos2)
    return pos2
end

function spin(angle::Real, pos::Vector)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.spin(angle, pos, pos2)
    return pos2
end

function ter2cel(jd_ut_high::Real,
    jd_ut_low::Real,
    delta_t::Real,
    method::Int,
    accuracy::Int,
    option::Int,
    xp::Real,
    yp::Real,
    vec::Vector)
    vec2 = zeros(Cdouble, 3)
    LibNOVAS.ter2cel(jd_ut_high, jd_ut_low, delta_t, method, accuracy, option, xp, yp, vec, vec)
    return vec2
end

function refract(location::LibNOVAS.on_surface, ref_option::Int, zd_obs::Real)
    LibNOVAS.refract(location, ref_option, zd_obs)
end

function equ2hor(jd_ut1::Real,
    delta_t::Real,
    accuracy::Int,
    xp::Real,
    yp::Real,
    location::NOVAS.OnSurface,
    ra::Real,
    dec::Real,
    ref_option::Int)
    zd = Ref{Cdouble}(0.0)
    az = Ref{Cdouble}(0.0)
    rar = Ref{Cdouble}(0.0)
    decr = Ref{Cdouble}(0.0)
    LibNOVAS.equ2hor(jd_ut1, delta_t, accuracy, xp, yp, location, ra, dec, ref_option, zd, az, rar, decr)
    return zd[], az[], rar[], decr[]
end