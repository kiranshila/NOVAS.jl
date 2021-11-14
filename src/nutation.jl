using DelimitedFiles,DataFrames

"""
"""


"""
    read_iau2000a()

Reads the included iau2000a nutation model into a DataFrame and returns them as
`(planetary_nutations,lunisolar_nutations)`
"""
function read_iau2000a()
    # Read planetary nutation terms
    lunisolar_nutations = DataFrame(readdlm("data/tab5.3a.txt";comment_char='*',comments=true)[1:678,:],
                                    [:l,:l′,:F,:D,:Ω,:Period,:A,:A′,:B,:B′,:A″,:A‴,:B″,:B‴])
    planetary_nutations = DataFrame(readdlm("data/tab5.3b.txt";skipstart=5),
                                    [:term,:l,:l′,:F,:D,:Ω,
                                    :Me,:Ve,:E,:Ma,:J,:Sa,:U,:Ne,:pA,
                                    :Period,:A,:A″,:B,:B″,:amplitude])
    return (planetary_nutations,lunisolar_nutations)
end

"""
    iau2000a(jd_high,jd_low)

Compute the forced nutation of the non-rigid earth based on the IAU 2000A nutation model.

# Arguments
-`jd_high::Real`: High order part of the TT Julian date
-`jd_low::Real`: Low order part of the TT Julian date

# Return
`(dpsi,deps)` where
-`dpsi`: Nutation (luni-solar + planetary) in longitude, in radians
-`deps`: Nutation (luni-solar + planetary) in obliquity, in radians
"""
function iau2000a(jd_high::Real,jd_low::Real)
    # How to cache this data?
    planetary,lunisolar = read_iau2000a()
    # Interval between fundamental epoch J2000.0 and given date
    t = ((jd_high-T0) + jd_low)/36525.0
    # Compute fundamental arguments in radians

    # Lunisolar nutation
    # nals_t is :l,:l′,:F,:D,Ω from lunisolar

    # Longitude
    Δψ = sum((lunisolar[!,:A]  + lunisolar[!,:A′]*t)*sin(argument) + 
             (lunisolar[!,:A″] + lunisolar[!,:A‴]*t)*cos(argument))
    # Obliquity
    Δε = sum((lunisolar[!,:B]  + lunisolar[!,:B′]*t)*sin(argument) + 
             (lunisolar[!,:B″] + lunisolar[!,:B‴]*t)*cos(argument))


end

export read_iau2000a,iau2000a
