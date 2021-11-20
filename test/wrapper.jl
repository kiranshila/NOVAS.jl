using NOVAS_jll

function fund_args(t::Real)
    a = zeros(Cdouble, 5)
    ccall((:fund_args, libnovas),
        Cvoid,
        (Cdouble, Ref{Cdouble}),
        t, a)
    return a
end

function norm_ang(angle::Real)
    ccall((:norm_ang, libnovas),
        Cdouble,
        (Cdouble,),
        angle)
end

function ee_ct(jd_high::Real, jd_low::Real, accuracy::Int)
    ccall((:ee_ct, libnovas),
        Cdouble,
        (Cdouble, Cdouble, Cshort),
        jd_high, jd_low, accuracy)
end

function iau2000a(jd_high::Real, jd_low::Real)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:iau2000a, libnovas),
        Cvoid,
        (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
        jd_high, jd_low, dpsi, deps)
    return dpsi[], deps[]
end

function nu2000k(jd_high::Real, jd_low::Real)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:nu2000k, libnovas),
        Cvoid,
        (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
        jd_high, jd_low, dpsi, deps)
    return dpsi[], deps[]
end

# Accuracy flags, 0 - Full, 1 - Reduced

function nutation_angles(t::Real, accuracy::Int)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:nutation_angles, libnovas),
        Cvoid,
        (Cdouble, Cshort, Ref{Cdouble}, Ref{Cdouble}),
        t, accuracy, dpsi, deps)
    return dpsi[], deps[]
end

function e_tilt(jd_tdb::Real, accuracy::Int)
    mobl = Ref{Cdouble}(0.0)
    tobl = Ref{Cdouble}(0.0)
    ee = Ref{Cdouble}(0.0)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:e_tilt, libnovas),
        Cvoid,
        (Cdouble, Cshort, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        jd_tdb, accuracy, mobl, tobl, ee, dpsi, deps)
    return mobl[], tobl[], ee[], dpsi[], deps[]
end

function mean_obliq(jd_tdb::Real)
    ccall((:mean_obliq, libnovas),
        Cdouble,
        (Cdouble,),
        jd_tdb)
end

function frame_tie(pos::AbstractVector,direction::Int)
    pos2 = zeros(Cdouble,3)
    ccall((:frame_tie, libnovas),
        Cvoid,
        (Ref{Cdouble},Cshort, Ref{Cdouble}),
        pos, direction, pos2)
    return pos2
end

function nutation(jd_tdb::Real,direction::Int,accuracy::Int,pos::AbstractVector)
    pos2 = zeros(Cdouble,3)
    ccall((:nutation, libnovas),
        Cvoid,
        (Cdouble, Cshort, Cshort,Ref{Cdouble}, Ref{Cdouble}),
        jd_tdb, direction, accuracy,pos,pos2)
    return pos2
end

function precession(jd_tdb1::Real,pos::AbstractVector,jd_tdb2::Real)
    pos2 = zeros(Cdouble,3)
    ccall((:precession, libnovas),
        Cshort,
        (Cdouble, Ref{Cdouble}, Cdouble, Ref{Cdouble}),
        jd_tdb1, pos, jd_tdb2,pos2)
    return pos2
end

function ira_equinox(jd_tdb::Real,equinox::Int,accuracy::Int)
    ccall((:ira_equinox,libnovas),
        Cdouble,
        (Cdouble,Cshort,Cshort),
        jd_tdb,equinox,accuracy)
end

function cio_location(jd_tdb::Real,accuracy::Int)
    ra_cio = Ref{Cdouble}(0.0)
    ref_sys = Ref{Cdouble}(0.0)
    ccall((:cio_location,libnovas),
        Cshort,
        (Cdouble,Cshort,Ref{Cdouble},Ref{Cdouble}),
        jd_tdb,accuracy,ra_cio,ref_sys)
    return ra_cio[]
end

function cio_basis(jd_tdb::Real,ra_cio::Real,ref_sys::Int,accuracy::Int)
    x = zeros(Cdouble,3)
    y = zeros(Cdouble,3)
    z = zeros(Cdouble,3)
    ccall((:cio_basis,libnovas),
        Cshort,
        (Cdouble,Cdouble,Cshort,Cshort,Ref{Cdouble},Ref{Cdouble},Ref{Cdouble}),
        jd_tdb,ra_cio,ref_sys,accuracy,x,y,z)
    return x,y,z
end

function era(jd_high::Real,jd_low::Real)
    ccall((:era,libnovas),
        Cdouble,
        (Cdouble,Cdouble),
        jd_high,jd_low)
end

function tdb2tt(jdb_jd::Real)
    tt_jd = Ref{Cdouble}(0.0)
    secdiff = Ref{Cdouble}(0.0)
    ccall((:tdb2tt,libnovas),
        Cvoid,
        (Cdouble,Ref{Cdouble},Ref{Cdouble}),
        jdb_jd,tt_jd,secdiff)
    return tt_jd[],secdiff[]
end