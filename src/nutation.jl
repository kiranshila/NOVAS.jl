using DelimitedFiles

include("utils.jl")

"""
    fund_args(t)

Computes the fundamental arguments (mean elements) of the Sun and Moon.

# Arguments
-`t::Real`: TDB time in Julian centuries since J2000.0

# Returns
`[l,l′,F,D,Ω]` where
-`l`: Mean anomaly of the Moon
-`l′`: Mean anomaly of the Sun
-`F`: mean argument of the latitude of the Moon
-`D`: Mean elongation of the Moon from the Sun
-`Ω` Mean longitude of the Moon's ascending node
"""
function fund_args(t::Real)
    l  = 485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470))))
    l′ = 1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149))))
    F  = 335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417))))
    D  = 1072260.70369 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169))))
    Ω  = 450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939))))
    return mod.([l, l′, F, D, Ω], ASEC360) .* ASEC2RAD
end

"""
    read_iau2000a()

Reads the included iau2000a nutation model into a Matrix and returns them as
`(planetary_nutations,lunisolar_nutations)`. This is meant to be a high-level interface for
working with the iau2000a data.
"""
function read_iau2000a()
    # Read planetary nutation terms
    ls_data = readdlm("data/tab5.3a.txt";comment_char='*',comments=true)[1:678,:]
    pl_data = readdlm("data/tab5.3b.txt";skipstart=5)
    return (pl_data, ls_data)
end

function iau2000a_statics()
    planetary, lunisolar = read_iau2000a()
    nals = lunisolar[:,1:5]                     # [:l,:l′,:F,:D,:Ω]
    cls = lunisolar[:,[7,8,11,9,10,13]] .* 1e4  # [:A, :A′, :A″, :B, :B′, :B″]
    napl = planetary[:,2:15]                    # [:l, :l′, :F, :D, :Ω, :Me, :Ve, :E, :Ma, :J, :Sa, :U, :Ne, :pA]
    cpl = planetary[:,17:20] .* 1e4             # [:A, :A″, :B, :B″]
    return nals, cls, napl, cpl
end

# Keep these around
const nals, cls, napl, cpl = iau2000a_statics()

"""
    iau2000a(jd)

Compute the forced nutation of the non-rigid earth based on the IAU 2000A nutation model.

# Arguments
-`jd_high::Real`: High order part of the TT Julian date
-`jd_low::Real=0.0`: Low order part of the TT Julian date

# Return
`(dpsi,deps)` where
-`dpsi`: Nutation (luni-solar + planetary) in longitude, in radians
-`deps`: Nutation (luni-solar + planetary) in obliquity, in radians
"""
function iau2000a(jd_high::Real, jd_low::Real=0.0)
    # Interval between fundamental epoch J2000.0 and given date
    t = ((jd_high - T0) + jd_low) / 36525.0

    # Compute fundamental arguments in radians
    a = fund_args(t)

    Δϕ_ls = 0.0
    Δε_ls = 0.0

    @inbounds @simd for i in 1:size(nals)[1]
        arg = nals[i,1] * a[1]  +
              nals[i,2] * a[2]  +
              nals[i,3] * a[3]  +
              nals[i,4] * a[4]  +
              nals[i,5] * a[5] 
        sarg, carg = sincos(arg)
        Δϕ_ls += (cls[i,1] + cls[i,2] * t) * sarg + cls[i,3] * carg;
        Δε_ls += (cls[i,4] + cls[i,5] * t) * carg + cls[i,6] * sarg;
    end

    # Mean anomaly of the Moon.
    al = mod2pi(2.35555598 + 8328.6914269554 * t)
    # Mean anomaly of the Sun.
    alsu = mod2pi(6.24006013 + 628.301955 * t)
    #  Mean argument of the latitude of the Moon.
    af = mod2pi(1.627905234 + 8433.466158131 * t)
    #  Mean elongation of the Moon from the Sun.
    ad = mod2pi(5.198466741 + 7771.3771468121 * t)
    # Mean longitude of the ascending node of the Moon.
    aom = mod2pi(2.18243920 - 33.757045 * t)
    # General accumulated precession in longitude.
    apa = (0.02438175 + 0.00000538691 * t) * t
    # Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
    alme = mod2pi(4.402608842 + 2608.7903141574 * t)
    alve = mod2pi(3.176146697 + 1021.3285546211 * t)
    alea = mod2pi(1.753470314 +  628.3075849991 * t)
    alma = mod2pi(6.203480913 +  334.0612426700 * t)
    alju = mod2pi(0.599546497 +   52.9690962641 * t)
    alsa = mod2pi(0.874016757 +   21.3299104960 * t)
    alur = mod2pi(5.481293871 +    7.4781598567 * t)
    alne = mod2pi(5.321159000 +    3.8127774000 * t)

    Δϕ_pl = 0.0
    Δε_pl = 0.0

    @inbounds @simd for i ∈ 1:size(napl)[1]
        arg = napl[i,1] * al    +
            napl[i,2] * alsu  +
            napl[i,3] * af    +
            napl[i,4] * ad    +
            napl[i,5] * aom   +
            napl[i,6] * alme  +
            napl[i,7] * alve  +
            napl[i,8] * alea  +
            napl[i,9] * alma  +
            napl[i,10] * alju  +
            napl[i,11] * alsa  +
            napl[i,12] * alur  +
            napl[i,13] * alne  +
            napl[i,14] * apa
        sarg, carg = sincos(arg)
        Δϕ_pl += cpl[i,1] * sarg + cpl[i,2] * carg;
        Δε_pl += cpl[i,3] * sarg + cpl[i,4] * carg;
    end

    # Scale Results
    scaling_factor = 1.0e-7 * ASEC2RAD
    
    # Return result!
    return (Δϕ_ls + Δϕ_pl, Δε_ls + Δε_pl) .* scaling_factor
end

export read_iau2000a,iau2000a,fund_args
