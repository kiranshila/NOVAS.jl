module NOVAS

include("constants.jl")
using LinearAlgebra

# Argument Enums
@enum Accuracy full reduced
@enum SiderealMethod equinox CIO
@enum GSTType mean_gst apparent_gst
@enum EquinoxType mean_equinox true_equinox
@enum NutationDirection mean2true true2mean
@enum FrameTieDirection dynamic2icrs icrs2dynamic

"""
    e_tilt(jd_tdb)

Computes quantities related to the orientation of the Earth's rotation axis at Julian date `jd_tdb`

# Arguments
- `jd_tdb::Real`: TDB Julian Date
- `accuracy::Accuracy=full`: Sets the accuracy level of `full` or `reduced`

# Returns
`(mobl,tobl,ee,dpsi,deps)` where
- `mobl`: Mean obliquity of the ecliptic in degrees
- `tobl`: True obliquity of the ecliptic in degrees
- `ee`: Equation of the equinoxes in seconds of time
- `dpsi`: Nutation in longitude on arcseconds
- `deps`: Nutation in obliquity in arcseconds
"""
# FIXME
function e_tilt(jd_tdb::Real;accuracy::Accuracy=full)
    # Compute time in Julian centuries from epoch J2000.0
    t = (jd_tdb - T0) / 36525.0
end

"""
    ira_equinox(jd_tdb)

Computes the intermediate right ascension of the equinox at the input Julian date using an
analytical expression for the accumulated precession in right ascension. For the true equinox,
the result is the equation of the origins.

# Arguments
- `jd_tdb::Real`: TDB Julian Date

# Optional Arguments
- `accuracy::Accuracy=full`: Either `full` or `reduced` accuracy
- `equinox::EquinoxType=mean_equinox`: Either `mean_equinox` or `true_equinox`
"""
function ira_equinox(jd_tdb::Real;accuracy::Accuracy=full,equinox::EquinoxType=mean_equinox)
    # NOTE: This function is memoized in novas.c
    # Compute time in Julian centuries
    t = (jd_tdb - T0 ) / 36525.0
    # For the true equinox, obtain the equation of the equinoxes in seconds,
    # which includes the 'complementary terms'
    if equinox == true_equinox
        _, _, eq_eq, _, _ = e_tilt(jd_tdb; accuracy=accuracy)
    elseif equinox == mean_equinox
        eq_eq = 0.0
    end
    # Compute precession in RA in arcseconds taken from the reference
    prec_ra = 0.014506 + ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t + 4612.156534) * t
    ra_eq = - (prec_ra / 15.0 + eq_eq) / 3600.0
    return ra_eq
end


"""
    nutation(jd_tdb)

Nutates equitorial rectangular coordinates from mean equator and equinox of epoch to true equator or
does the transformation in reverse given the flag `direction`.

# Arguments
- `jd_tdb::Real`: TDB Julian Date
- `pos::AbstractVector`: Position vector, geocentric equatorial rectangular coordinates, referred to mean equator and equinox of epoch

# Optional Arguments
- `accuracy::Accuracy=full`: Either `full` or `reduced` accuracy
- `direction::NutationDirection=mean2true`: Perform the transformation in either the `mean2true` or `true2mean` direction
"""
function nutation(jd_tdb::Real, pos::AbstractVector;accuracy::Accuracy=full,direction::NutationDirection=mean2true)
    # Call e_tilt to get the obliquity and nutation angles
    oblm, oblt, _, psi, _ = e_tilt(jd_tdb;accuracy=accuracy)
    # Construct rotation matrix
    sobm, cobm = sincosd(oblm)
    sobt, cobt = sincosd(oblt)
    spsi, cpsi = sincos(psi * ASEC2RAD)
    rot = [  cpsi            -spsi * cobm                -spsi * sobm         ;
           spsi * cobt   cpsi * cobm * cobt + sobm * sobt   cpsi * sobm * cobt - cobm * sobt;
           spsi * sobt   cpsi * cobm * sobt - sobm * cobt   cpsi * sobm * sobt + cobm * cobt]
    # Transform Position
    if direction == mean2true
        return rot * pos
    elseif direction == true2mean
        return rot' * pos
    end
end

"""
    precession(jd_tdb1,pos,jd_tdb2)

Precesses equatorial rectangular coordinates from one epoch to another. One of the two epochs must be J200.0.

# Arguments
- `jd_tdb1::Real`: TDB Julian Date of first epoch
- `jd_tdb2::Real`: TDB Julian Date of second epoch
- `pos::AbstractVector`: Position vector, geocentric equatorial rectangular coordinates, referred to mean dynamical equator and equinox of first epoch
"""
function precession(jd_tdb1::Real, pos::AbstractVector, jd_tdb2::Real)
    # NOTE: this function is memoized in novas.c
    @assert (jd_tdb1 == T0) || (jd_tdb2 == T0) "One of jd_tdb1, jd_tdb2 must be T0"
    eps0 = 84381.406
    # t is time in TDB centuries between the two epochs
    t = (jd_tdb2 - jd_tdb1) / 36525.0
    if jd_tdb2 == T0
        t *= -1
    end

    psia   = ((((-    0.0000000951  * t
                   +    0.000132851 ) * t
                   -    0.00114045  ) * t
                   -    1.0790069   ) * t
                   + 5038.481507    ) * t
    omegaa = ((((+    0.0000003337  * t
                   -    0.000000467 ) * t
                   -    0.00772503  ) * t
                   +    0.0512623   ) * t
                   -    0.025754    ) * t + eps0
    chia   = ((((-    0.0000000560  * t
                   +    0.000170663 ) * t
                   -    0.00121197  ) * t
                   -    2.3814292   ) * t
                   +   10.556403    ) * t

    eps0 = eps0 * ASEC2RAD
    psia = psia * ASEC2RAD
    omegaa = omegaa * ASEC2RAD
    chia = chia * ASEC2RAD

    sa, ca = sincos(eps0)
    sb, cb = sincos(-psia)
    sc, cc = sincos(-omegaa)
    sd, cd = sincos(chia)

    # Construct precession rotation matrix
    rot = [ cd * cb - sb * sd * cc   cd * sb * ca + sd * cc * cb * ca - sa * sd * sc   cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
           -sd * cb - sb * cd * cc  -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc  -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
                 sb * sc             -sc * cb * ca - sa * cc                 -sc * cb * sa + cc * ca]

    # Perform rotation from epoch to J2000.0
    if jd_tdb2 == T0
        return rot' * pos
    else
        return rot * pos
    end
end

"""
    frame_tie(pos;direction)

Transforms a position vector from the dynamical reference system to ICRS, or vice versa.

# Argumenets
- `pos::AbstractVector`: Position vector, equitorial rectangular coordinates
- `direction::FrameTieDirection=dynamic2icrs`
"""
function frame_tie(pos::AbstractVector, direction::FrameTieDirection=dynamic2icrs)
    if direction == dynamic2icrs
        return frame_tie_rot * pos
    else
        return frame_tie_rot' * pos
    end
end

"""
   cio_location(jd_tdb)

Computes the location of the celestial intermediate origin (CIO) for a given Julian date, as a
right ascension with respect to the true equinox of date.

# Arguments
- `jd_tdb::Real`: TDB Julian Date

# Optional Arguments
- `accuracy::Accuracy=full`: Either `full` or `reduced` accuracy
"""
function cio_location(jd_tdb::Real;accuracy::Accuracy=full)
    # NOTE: This function is memoized in novas.c
    # NOTE: We are going to calculate the values instead of using the external dep for the time being
    ra_cio = -ira_equinox(jd_tdb, accuracy=accuracy)
    return ra_cio
end

"""
    cio_basis(jd_tdb,ra_cio)

Computes the orthonomal basis vectors, with respect to the GCRS of the celestial intermediate system defined
by the celestial intermediate pole and origin.

# Arguments
- `jd_tdb::Real`: TDB Julian Date

# Optional Arguments
- `accuracy::Accuracy=full`: Either `full` or `reduced` accuracy
"""
# FIXME
function cio_basis(jd_tdb::Real, ra_cio::Real;accuracy::Accuracy=full)
    # NOTE: This function is memoized in novas.c

    # Compute unit vector z towards celestial pole
    z0 = [0.0, 0.0, 1.0]
    w1 = nutation(jd_tdb, -1, accuracy, z0)
    w2 = precession(jd_tdb, w1, T0)
    zz = frame_tie(w2, 01)
end

"""
    era(jd_high)

This function returns the value of the earth rotation angle (theta in degrees) for a given UT1 Julian date.

# Arguments
- `jd_high::Real`: High-order part of UT1 Julian date
- `jd_low::Real`: Low-order part of UT1 Julian date
"""
function era(jd_high::Real, jd_low::Real)
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
    return (tt_jd, secdiff)
end

"""
    sidereal_time(jd_ut1,0,delta_t)

Computes the Greenwich sidereal time, either mean or apparent, at Julian date `jd_high` + `jd_low` in hours

# Arguments
- `jd_high::Real`: high-order part of UT1 Julian date
- `jd_low::Real`: low-order part of UT1 Julian date
- `delta_t::Real`: Difference TT-UT1 at `jd_high` + `jd_low` in seconds

# Optional arguments
- `gst_type::GSTType=mean_gst`: Return results as mean (`mean_gst`) or apparent (`apparent_gst`) time
- `method::SiderealMethod=equinox`: Computation method, CIO-based (`CIO`) or equinox-based (`equinox`)
- `accuracy::Accuracy=full`: Either `full` or `reduced` accuracy
"""
function sidereal_time(jd_high::Real,
                       jd_low::Real,
                       delta_t::Real;
                       gst_type::GSTType=mean_gst,
                       method::SiderealMethod=equinox,
                       accuracy::Accuracy=full)
    # NOTE: This function is memoized in novas.c

    # Time argument for precession and nutation components of sidereal time is TDB. First approximation
    # is TDB = TT, then refine
    jd_ut = jd_high + jd_low
    jd_tt = jd_ut + (delta_t / 86400.0)
    jd_tdb = jd_tt
    _, secdiff = tdb2tt(jd_tdb)
    jd_tdb = jd_tt + (secdiff / 86400.0)
    t = (jd_tdb - T0) / 36525.0

    # Compute the Earth Rotation Angle. Time argument is UT1.
    theta = era(jd_high, jd_low)

    # Compute the equation of the equinoxes if needed, depending upon the input values of `gst_type` and `method`.
    # If not needed, set to zero
    # FIXME
    # if ((gst_type == :mean)     && (method == :CIO)) ||
    #    ((gst_type == :apparent) && (method == :equinox))
    #     _,_,ee,_,_ = e_tilt(jd_tdb,accuracy)
    #     eqeq = ee * 15
    # else
    #     eqeq = 0
    # end
    eqeq = 0

    if method == CIO
        # Obtain basis vectors, in the GCRS of the celestial intermediate system
        ra_cio = cio_location(jd_tdb, accuracy=accuracy)
    elseif method == equinox
        # Precession-in-RA terms in mean sidereal time with coefficients in arcseconds
        st = eqeq + 0.014506 + ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t + 4612.156534) * t
        gst = mod((st / 3600.0 + theta), 360.0) / 15.0;
        if gst < 0.0
            gst += 24.0
        end
    end
    return gst
end

# Enum Exports
export GSTType,SiderealMethod,Accuracy,EquinoxType, FrameTieDirection
# Function exports
export sidereal_time,tdb2tt,era,e_tilt,ira_equinox,nutation,precession,frame_tie,cio_location,cio_basis

end
