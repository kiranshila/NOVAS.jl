using NOVAS
import NOVAS.LibNOVAS

function fund_args(t)
    a = zeros(Cdouble, 5)
    LibNOVAS.fund_args(t, a)
    return a
end

function norm_ang(angle)
    return LibNOVAS.norm_ang(angle)
end

function ee_ct(jd_high, jd_low; accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    return LibNOVAS.ee_ct(jd_high, jd_low, acc)
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

function nutation_angles(t; accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.nutation_angles(t, acc, dpsi, deps)
    return dpsi[], deps[]
end

function e_tilt(jd_tdb; accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    mobl = Ref{Cdouble}(0.0)
    tobl = Ref{Cdouble}(0.0)
    ee = Ref{Cdouble}(0.0)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    LibNOVAS.e_tilt(jd_tdb, acc, mobl, tobl, ee, dpsi, deps)
    return mobl[], tobl[], ee[], dpsi[], deps[]
end

function mean_obliq(jd_tdb)
    return LibNOVAS.mean_obliq(jd_tdb)
end

function frame_tie(pos, direction)
    if direction == :dynamic2icrs
        dir = -1
    else
        dir = 0
    end
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.frame_tie(pos, dir, pos2)
    return pos2
end

function nutation(jd_tdb, pos; direction, accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end

    if direction == :mean2true
        dir = 0
    elseif direction == :true2mean
        dir = 1
    end

    pos2 = zeros(Cdouble, 3)
    LibNOVAS.nutation(jd_tdb, dir, acc, pos, pos2)
    return pos2
end

function precession(jd_tdb1, pos, jd_tdb2)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.precession(jd_tdb1, pos, jd_tdb2, pos2)
    return pos2
end

function ira_equinox(jd_tdb; equinox, accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    if equinox == :mean
        equ = 0
    elseif equinox == :true
        equ = 1
    end
    return LibNOVAS.ira_equinox(jd_tdb, equ, acc)
end

function cio_location(jd_tdb; accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    ra_cio = Ref{Cdouble}(0.0)
    ref_sys = Ref{Cshort}(0)
    LibNOVAS.cio_location(jd_tdb, acc, ra_cio, ref_sys)
    return ra_cio[]
end

function cio_basis(jd_tdb, ra_cio; accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    x = zeros(Cdouble, 3)
    y = zeros(Cdouble, 3)
    z = zeros(Cdouble, 3)
    LibNOVAS.cio_basis(jd_tdb, ra_cio, 2, acc, x, y, z)
    return x, y, z
end

function era(jd_high, jd_low)
    return LibNOVAS.era(jd_high, jd_low)
end

function tdb2tt(jdb_jd)
    tt_jd = Ref{Cdouble}(0.0)
    secdiff = Ref{Cdouble}(0.0)
    LibNOVAS.tdb2tt(jdb_jd, tt_jd, secdiff)
    return tt_jd[], secdiff[]
end

function sidereal_time(jd_high, jd_low, delta_t; gst_type, method, accuracy)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end

    if gst_type == :mean
        gstt = 0
    elseif gst_type == :apparent
        gstt = 1
    end

    if method == :CIO
        meth = 0
    elseif method == :equinox
        meth = 1
    end

    gst = Ref{Cdouble}(0.0)
    LibNOVAS.sidereal_time(jd_high, jd_low, delta_t, gstt, meth, acc, gst)
    return gst[]
end

function wobble(tjd, xp, yp, pos; direction)
    if direction == :itrs2terr
        dir = 0
    elseif direction == :terr2itrs
        dir = 1
    end
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.wobble(tjd, dir, xp, yp, pos, pos2)
    return pos2
end

function spin(angle, pos)
    pos2 = zeros(Cdouble, 3)
    LibNOVAS.spin(angle, pos, pos2)
    return pos2
end

function ter2cel(jd_ut_high, jd_ut_low, delta_t, pos; method, accuracy, option, xp, yp)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    if method == :CIO
        meth = 0
    elseif method == :equinox
        meth = 1
    end
    if option == :GCRS
        opt = 0
    elseif option == :equinox
        opt = 1
    end
    vec2 = zeros(Cdouble, 3)
    LibNOVAS.ter2cel(jd_ut_high, jd_ut_low, delta_t, meth, acc, opt, xp, yp, pos, vec2)
    return vec2
end

function refract(location, zd_obs; ref_option)
    if ref_option == :standard
        ro = 1
    elseif ref_option == :location
        ro = 2
    end
    loc = convert(NOVAS.LibNOVAS.on_surface, location)
    return LibNOVAS.refract(Ref{NOVAS.LibNOVAS.on_surface}(loc), ro, zd_obs)
end

function equ2hor(jd_ut1, delta_t, ra, dec, location; accuracy, xp, yp, ref_option)
    zd = Ref{Cdouble}(0.0)
    az = Ref{Cdouble}(0.0)
    rar = Ref{Cdouble}(0.0)
    decr = Ref{Cdouble}(0.0)
    if accuracy == :full
        acc = 0
    elseif accuracy == :reduced
        acc = 1
    end
    if ref_option == :none
        ro = 0
    elseif ref_option == :standard
        ro = 1
    elseif ref_option == :location
        ro = 2
    end
    loc = Ref{NOVAS.LibNOVAS.on_surface}(convert(NOVAS.LibNOVAS.on_surface, location))
    LibNOVAS.equ2hor(jd_ut1, delta_t, acc, xp, yp, loc, ra, dec, ro, zd, az, rar, decr)
    return zd[], az[], rar[], decr[]
end

import Base.convert

function Base.convert(::Type{NOVAS.LibNOVAS.on_surface}, obj::OnSurface{Float64})
    # There's gotta be a better way to do this
    return NOVAS.LibNOVAS.on_surface(obj.latitude, obj.longitude, obj.height,
                                     obj.temperature, obj.pressure)
end
