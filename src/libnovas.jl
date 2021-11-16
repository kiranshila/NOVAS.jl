# This defines the C wrapping interface
using Libdl

function novas_iau2000a(jd_high::Real,jd_low::Real)
    dpsi = Ref{Cdouble}(0)
    deps = Ref{Cdouble}(0)
    iau2000a_sym = Libdl.dlsym(libnovas, :iau2000a)
    ccall(iau2000a_sym ,
        Cvoid,
        (Cdouble,Cdouble,Ref{Cdouble},Ref{Cdouble}),
        jd_high,jd_low,dpsi,deps)
    return dpsi[],deps[]
end
