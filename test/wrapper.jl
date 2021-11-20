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