using DelimitedFiles

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

"""
    read_nu2000k()

Reads the included nu2000k nutation model into a Matrix and returns them as
`(planetary_nutations,lunisolar_nutations)`. This is meant to be a high-level interface for
working with the nu2000k data.
"""
function read_nu2000k()
    napl = readdlm("data/nu2000k_napl.csv",',')
    nals = readdlm("data/nu2000k_nals.csv",',')
    cpl = readdlm("data/nu2000k_cpl.csv",',')
    cls = readdlm("data/nu2000k_cls.csv",',')
    return nals, cls, napl, cpl
end

# Keep these around
const nals_a, cls_a, napl_a, cpl_a = iau2000a_statics()
const nals_k, cls_k, napl_k, cpl_k = read_nu2000k()

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

    @inbounds @simd for i in 1:size(nals_a)[1]
        arg = nals_a[i,1] * a[1]  +
              nals_a[i,2] * a[2]  +
              nals_a[i,3] * a[3]  +
              nals_a[i,4] * a[4]  +
              nals_a[i,5] * a[5] 
        sarg, carg = sincos(arg)
        Δϕ_ls += (cls_a[i,1] + cls_a[i,2] * t) * sarg + cls_a[i,3] * carg
        Δε_ls += (cls_a[i,4] + cls_a[i,5] * t) * carg + cls_a[i,6] * sarg
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

    @inbounds @simd for i ∈ 1:size(napl_a)[1]
        arg = napl_a[i,1]  * al    +
              napl_a[i,2]  * alsu  +
              napl_a[i,3]  * af    +
              napl_a[i,4]  * ad    +
              napl_a[i,5]  * aom   +
              napl_a[i,6]  * alme  +
              napl_a[i,7]  * alve  +
              napl_a[i,8]  * alea  +
              napl_a[i,9]  * alma  +
              napl_a[i,10] * alju  +
              napl_a[i,11] * alsa  +
              napl_a[i,12] * alur  +
              napl_a[i,13] * alne  +
              napl_a[i,14] * apa
        sarg, carg = sincos(arg)
        Δϕ_pl += cpl_a[i,1] * sarg + cpl_a[i,2] * carg
        Δε_pl += cpl_a[i,3] * sarg + cpl_a[i,4] * carg
    end

    # Scale Results
    scaling_factor = 1.0e-7 * ASEC2RAD
    
    # Return result!
    return (Δϕ_ls + Δϕ_pl, Δε_ls + Δε_pl) .* scaling_factor
end


"""
    nu2000k(jd)

Compute the forced nutation of the non-rigid earth based on the NU2000K nutation model.

# Arguments
-`jd_high::Real`: High order part of the TT Julian date
-`jd_low::Real=0.0`: Low order part of the TT Julian date

# Return
`(dpsi,deps)` where
-`dpsi`: Nutation (luni-solar + planetary) in longitude, in radians
-`deps`: Nutation (luni-solar + planetary) in obliquity, in radians
"""
function nu2000k(jd_high::Real, jd_low::Real=0.0)
    # Interval between fundamental epoch J2000.0 and given date
    t = ((jd_high - T0) + jd_low) / 36525.0

    # Compute fundamental arguments in radians
    a = fund_args(t)

    Δϕ_ls = 0.0
    Δε_ls = 0.0

    @inbounds @simd for i in 1:size(nals_k)[1]
        arg = nals_k[i,1] * a[1]  +
              nals_k[i,2] * a[2]  +
              nals_k[i,3] * a[3]  +
              nals_k[i,4] * a[4] + 
              nals_k[i,5] * a[5]
        sarg, carg = sincos(arg)
        Δϕ_ls += (cls_k[i,1] + cls_k[i,2] * t) * sarg + cls_k[i,3] * carg
        Δε_ls += (cls_k[i,4] + cls_k[i,5] * t) * carg + cls_k[i,6] * sarg
    end

    # Planetary longitudes, Mercury through Neptune, wrt mean dynamical
    # ecliptic and equinox of J2000, with high order terms omitted
    # (Simon et al. 1994, 5.8.1-5.8.8).
    alme = mod2pi(4.402608842461 + 2608.790314157421 * t)
    alve = mod2pi(3.176146696956 + 1021.328554621099 * t)
    alea = mod2pi(1.753470459496 +  628.307584999142 * t)
    alma = mod2pi(6.203476112911 +  334.061242669982 * t)
    alju = mod2pi(0.599547105074 +   52.969096264064 * t)
    alsa = mod2pi(0.874016284019 +   21.329910496032 * t)
    alur = mod2pi(5.481293871537 +    7.478159856729 * t)
    alne = mod2pi(5.311886286677 +    3.813303563778 * t)
    # General precession in longitude (Simon et al. 1994), equivalent to 5028.8200 arcsec/cy at J2000.
    apa = (0.024380407358 + 0.000005391235 * t) * t

    Δϕ_pl = 0.0
    Δε_pl = 0.0

    @inbounds @simd for i ∈ 1:size(napl_k)[1]
        arg = napl_k[i,1]  * a[1] +
              napl_k[i,2]  * a[2] +
              napl_k[i,3]  * a[3] +
              napl_k[i,4]  * a[4] +
              napl_k[i,5]  * a[5] +
              napl_k[i,6]  * alme +
              napl_k[i,7]  * alve +
              napl_k[i,8]  * alea +
              napl_k[i,9]  * alma +
              napl_k[i,10] * alju +
              napl_k[i,11] * alsa +
              napl_k[i,12] * alur +
              napl_k[i,13] * alne +
              napl_k[i,14] * apa
        sarg, carg = sincos(arg)
        Δϕ_pl += cpl_k[i,1] * sarg + cpl_k[i,2] * carg;
        Δε_pl += cpl_k[i,3] * sarg + cpl_k[i,4] * carg;
    end

    # Scale Results
    scaling_factor = 1.0e-7 * ASEC2RAD
    
    # Return result!
    return (Δϕ_ls + Δϕ_pl, Δε_ls + Δε_pl) .* scaling_factor
end

export read_iau2000a,read_nu2000k,iau2000a,nu2000k,fund_args