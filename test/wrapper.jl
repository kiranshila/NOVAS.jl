using NOVAS_jll

function fund_args(t::Float64)
    a = zeros(Cdouble,5)
    ccall((:fund_args,libnovas),
          Cvoid,
          (Cdouble,Ref{Cdouble}),
          t,a)
    return a
end

function norm_ang(angle::Float64)
    ccall((:norm_ang,libnovas),
          Cdouble,
          (Cdouble,),
          angle)
end

function ee_ct(jd_high::Float64,jd_low::Float64,accuracy::Int)
    ccall((:ee_ct,libnovas),
        Cdouble,
        (Cdouble,Cdouble,Cshort),
        jd_high,jd_low,accuracy)
end

function iau2000a(jd_high::Float64,jd_low::Float64)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:iau2000a,libnovas),
          Cvoid,
          (Cdouble,Cdouble,Ref{Cdouble},Ref{Cdouble}),
          jd_high,jd_low,dpsi,deps)
    return dpsi[],deps[]
end

function nu2000k(jd_high::Float64,jd_low::Float64)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:nu2000k,libnovas),
          Cvoid,
          (Cdouble,Cdouble,Ref{Cdouble},Ref{Cdouble}),
          jd_high,jd_low,dpsi,deps)
    return dpsi[],deps[]
end

# Accuracy flags, 0 - Full, 1 - Reduced

function nutation_angles(t::Float64,accuracy::Int)
    dpsi = Ref{Cdouble}(0.0)
    deps = Ref{Cdouble}(0.0)
    ccall((:nutation_angles,libnovas),
        Cvoid,
        (Cdouble,Cshort,Ref{Cdouble},Ref{Cdouble}),
        t,accuracy,dpsi,deps)
    return dpsi[],deps[]
end

function e_tilt(jd_tdb::Real,accuracy::Int)


end