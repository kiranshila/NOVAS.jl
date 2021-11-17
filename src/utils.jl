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
    l = 485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470))))
    l′ = 1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149))))
    F = 335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417))))
    D = 1072260.70369 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169))))
    Ω = 450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939))))
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
    ls_data = readdlm("data/tab5.3a.txt"; comment_char = '*', comments = true)[1:678, :]
    pl_data = readdlm("data/tab5.3b.txt"; skipstart = 5)
    return (pl_data, ls_data)
end

"""
    read_nu2000k()

Reads the included nu2000k nutation model into a Matrix and returns them as
`(planetary_nutations,lunisolar_nutations)`. This is meant to be a high-level interface for
working with the nu2000k data.
"""
function read_nu2000k()
    napl = readdlm("data/nu2000k_napl.csv", ',')
    nals = readdlm("data/nu2000k_nals.csv", ',')
    cpl = readdlm("data/nu2000k_cpl.csv", ',')
    cls = readdlm("data/nu2000k_cls.csv", ',')
    return nals, cls, napl, cpl
end

"""
    read_cterms()

Reads the included complementary term constants.
"""
function read_cterms()
    cterms = readdlm("data/cterms.txt")
    ke0 = cterms[1:33, 4:17]
    ke1 = cterms[34, 4:17]
    se0 = cterms[1:33, 2:3] .* 1e-6
    se1 = cterms[34, 2:3] .* 1e-6
    return ke0, ke1, se0, se1
end

# Keep these around
const ke0, ke1, se0, se1 = read_cterms()


"""
    norm_ang(θ)

Normalizes an angle into range 0 <= θ <= 2π
"""
function norm_ang(angle::Real)
    a = mod2pi(angle);
    if a < 0.0
          a += TWOPI
    end
    return a
end

"""
    ee_ct(jd)

Computes the "complementary terms" of the equation of the equinoxes.

# Arguments
- `jd_high::Real`: High-order part of UT1 Julian date
- `jd_low::Real=0.0`: Low-order part of UT1 Julian date

# Optional Arguments
- `accuracy::Symbol=:full`: Sets the accuracy level of `:full` or `:reduced`
"""
function ee_ct(jd_high::T, jd_low::T = 0.0; accuracy::Symbol = :full) where {T<:Real}
    t = ((jd_high - T0) + jd_low) / 36525.0

    fa = zeros(T, 14)

    if accuracy == :full
        # Fundamental arguments
        # Mean anomaly of the moon
        fa[1] = norm_ang((485868.249036 + (715923.2178 + (31.8792 + (0.051635 + (-0.00024470) * t) * t) * t) * t) * ASEC2RAD + mod(1325.0 * t, 1.0) * 2π)
        # Mean anomaly of the Sun
        fa[2] = norm_ang((1287104.793048 +
                         (1292581.0481 +
                         (     -0.5532 +
                         (     +0.000136 +
                         (     -0.00001149)
                         * t) * t) * t) * t) * ASEC2RAD
                         + mod(99.0*t, 1.0) * 2π)
        # Mean Longitude of the Moon minus Mean Longitude of the Ascending Node of the Moon.
        fa[3] = norm_ang(( 335779.526232 +
                         ( 295262.8478 +
                         (    -12.7512 +
                         (     -0.001037 +
                         (      0.00000417)
                         * t) * t) * t) * t) * ASEC2RAD
                         + mod(1342.0*t, 1.0) * 2π)
        # Mean Elongation of the Moon from the Sun.
        fa[4] = norm_ang((1072260.703692 +
                         (1105601.2090 +
                         (     -6.3706 +
                         (      0.006593 +
                         (     -0.00003169)
                         * t) * t) * t) * t) * ASEC2RAD
                         + mod(1236.0*t, 1.0) * 2π)
        # Mean Longitude of the Ascending Node of the Moon.
        fa[5] = norm_ang(( 450160.398036 +
                         (-482890.5431 +
                         (      7.4722 +
                         (      0.007702 +
                         (     -0.00005939)
                         * t) * t) * t) * t) * ASEC2RAD
                         + mod(-5.0*t, 1.0) * 2π)
        fa[6] = norm_ang(4.402608842 + 2608.7903141574 * t)
        fa[7] = norm_ang(3.176146697 + 1021.3285546211 * t)
        fa[8] = norm_ang(1.753470314 +  628.3075849991 * t)
        fa[9] = norm_ang(6.203480913 +  334.0612426700 * t)
        fa[10] = norm_ang(0.599546497 +   52.9690962641 * t)
        fa[11] = norm_ang(0.874016757 +   21.3299104960 * t)
        fa[12] = norm_ang(5.481293872 +    7.4781598567 * t)
        fa[13] = norm_ang(5.311886287 +    3.8133035638 * t)
        fa[14] =          (0.024381750 +    0.00000538691 * t) * t
        # Evaluate the complementary terms.
        s0 = 0.0
        s1 = 0.0
        @inbounds @simd for i ∈ 1:33
            a = 0.0
            for j ∈ 1:14
                a += ke0[i,j] * fa[j]
            end
            s0 += se0[i,1] * sin(a) + se0[i,2] * cos(a)
        end
        a = 0.0
        @inbounds @simd for j ∈ 1:14
            a += ke1[j] * fa[j];
        end
        s1 += se1[1] * sin(a) + se1[2] * cos(a)
        c_terms = s0 + s1 * t

    elseif accuracy == :reduced
        fa2 = fund_args(t)
        c_terms =
            2640.96e-6 * sin(fa2[5])
           +  63.52e-6 * sin(2.0 * fa2[5])
           +  11.75e-6 * sin(2.0 * fa2[3] - 2.0 * fa2[4] + 3.0 * fa2[5])
           +  11.21e-6 * sin(2.0 * fa2[3] - 2.0 * fa2[4] +       fa2[5])
           -   4.55e-6 * sin(2.0 * fa2[3] - 2.0 * fa2[4] + 2.0 * fa2[5])
           +   2.02e-6 * sin(2.0 * fa2[3]                + 3.0 * fa2[5])
           +   1.98e-6 * sin(2.0 * fa2[3]                +       fa2[5])
           -   1.72e-6 * sin(3.0 * fa2[5])
           -   0.87e-6 * t * sin(fa2[5]);
    end
    return c_terms * ASEC2RAD
end

export ee_ct