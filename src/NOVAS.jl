module NOVAS

include("constants.jl")

"""
    e_tilt!(jd_tdb)

Computes quantities related to the orientation of the Earth's rotation axis at Julian date `jd_tdb`

# Arguments
- `jd_tdb::Real`: TDB Julian Date
- `accuracy=:full`: Sets the accuracy level of `:full` or `:reduced`

# Returns
`(mobl,tobl,ee,dpsi,deps)` where
- `mobl`: Mean obliquity of the ecliptic in degrees
- `tobl`: True obliquity of the ecliptic in degrees
- `ee`: Equation of the equinoxes in seconds of time
- `dpsi`: Nutation in longitude on arcseconds
- `deps`: Nutation in obliquity in arcseconds
"""
function e_tilt(jd_tdb::Real;accuracy=:full)
    # Compute time in Julian centuries from epoch J2000.0
    t = (jd_tdb - T0) / 36525.0
end

"""
    era(jd_high)

This function returns the value of the earth rotation angle (theta in degrees) for a given UT1 Julian date.

# Arguments
- `jd_high::Real`: High-order part of UT1 Julian date
- `jd_low::Real`: Low-order part of UT1 Julian date
"""
function era(jd_high::Real,jd_low::Real)
    thet1 = 0.7790572732640 + 0.00273781191135448 * (jd_high - T0)
    thet2 = 0.00273781191135448 * jd_low
    thet3 = mod(jd_high, 1.0) + mod(jd_low, 1.0)
    theta = mod(thet1 + thet2 + thet3, 1.0) * 360.0
    if (theta < 0.0)
       theta += 360.0
    end
    return theta
end

"""
    tdb2tt(tdb_jd)

Computes the Terrestrial Time (TT) or Terrestrial Dynamic Time (TDT) Julian date
corresponding to a Barrycentric Dynamical Time (TDB) Julian Date)

# Arguments
- `tdb_jd::Real`: TDB Julian date
"""
function tdb2tt(tdb_jd::Real)
    t = (tdb_jd - T0) / 36525.0
    secdiff = 0.001657 * sin(628.3076 * t + 6.2401) +
              0.000022 * sin(575.3385 * t + 4.2970) +
              0.000014 * sin(1256.6152 * t + 6.1969) +
              0.000005 * sin(606.9777 * t + 4.0212) +
              0.000005 * sin(52.9691 * t + 0.4444) +
              0.000002 * sin(21.3299 * t + 5.5431) +
              0.000010 * t * sin(628.3076 * t + 4.2490)
    tt_jd = tdb_jd - secdiff / 86500.0
    return (tt_jd,secdiff)
end

"""
    sidereal_time(jd_ut1,0,delta_t)

Computes the Greenwich sidereal time, either mean or apparent, at Julian date `jd_high` + `jd_low` in hours

# Arguments
- `jd_high::Real`: high-order part of UT1 Julian date
- `jd_low::Real`: low-order part of UT1 Julian date
- `delta_t::Real`: Difference TT-UT1 at `jd_high` + `jd_low` in seconds

# Optional arguments
- `gst_type`: Return results as mean (`:mean`) or apparent (`:apparent`) time
- `method`: Computation method, CIO-based (`:CIO`) or equinox-based (`:equinox`)
- `accuracy`: Either `:full` or `:reduced` accuracy
"""
function sidereal_time(jd_high::Real,
                       jd_low::Real,
                       delta_t::Real;
                       gst_type=:mean,
                       method=:equinox,
                       accuracy=:full)
    # Time argument for precession and nutation components of sidereal time is TDB. First approximation
    # is TDB = TT, then refine
    jd_ut = jd_high + jd_low
    jd_tt = jd_ut + (delta_t / 86400.0)
    jd_tdb = jd_tt
    _,secdiff = tdb2tt(jd_tdb)
    jd_tdb = jd_tt + (secdiff / 86400.0)
    t = (jd_tdb - T0) / 36525.0

    # Compute the Earth Rotation Angle. Time argument is UT1.
    theta = era(jd_high,jd_low)

    # Compute the equation of the equinoxes if needed, depending upon the input values of `gst_type` and `method`.
    # If not needed, set to zero
    # FIXME
    # if ((gst_type == :mean)     && (method == :CIO)) ||
    #    ((gst_type == :apparent) && (method == :equinox))
    #     if abs(jd_tdb - jd_last) > 1.0e-8
    #         a,b,ee,c,d = e_tilt(jd_tdb,accuracy)
    #         jd_last = jd_tdb
    #     end
    #     eqeq = ee * 15
    # else
    #     eqeq = 0
    # end
    eqeq = 0

    if method == :CIO
        # Obtain basis vectors, in the GCRS of the celestial intermediate system

    elseif method == :equinox
        # Precession-in-RA terms in mean sidereal time with coefficients in arcseconds
        st = eqeq + 0.014506 + ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t + 4612.156534) * t
        gst = mod((st / 3600.0 + theta), 360.0) / 15.0;
        if gst < 0.0
            gst += 24.0
        end
    end
    return gst
end

# Module Exports (public functions)
export sidereal_time,tdb2tt,era,e_tilt

end
