export nutation_angles, mean_obliq, frame_tie, nutation, precession, cel_pole, e_tilt,
       ira_equinox, cio_location, cio_basis, era, tdb2tt, sidereal_time, wobble, spin,
       ter2cel, OnSurface, refract, equ2hor

"""
    nutation_angles(t)

Computes the values for nutation in longitude and nutation in obliquity for a given TDB julian date.

# Aruments
- `t::Real`: TDB time in Julian centuries since J2000.0

# Optional Arguments
- `accuracy::Symbol=:full`: Sets the accuracy level of `:full` or `:reduced`

# Returns
`(dpsi,deps)` where
- `dpsi`: Nutation in longitude in arcseconds
- `deps`: Nutation in obliquity in arcseconds
"""
function nutation_angles(t::Real; accuracy::Symbol=:full)
    @assert accuracy ∈ Set([:full, :reduced])
    t1 = t * 36525.0
    # High accuracy mode uses IAU 200A
    if accuracy == :full
        dpsi, deps = iau2000a(T0, t1)
    elseif accuracy == :reduced
        # Low accuracy mode uses the specially truncated version of IAU 2000A, called NU200K
        dpsi, deps = nu2000k(T0, t1)
    end
    # Convert output to arcseconds
    dpsi /= ASEC2RAD
    deps /= ASEC2RAD
    # Return
    return (dpsi, deps)
end

"""
    mean_obliq(jd_tdb)

Computes the mean obliquity of the ecliptic
"""
function mean_obliq(jd_tdb::Real)
    # Julian centuries since J2000.0 epoch
    t = (jd_tdb - T0) / 36525.0
    ε = ((((-0.0000000434 * t - 0.000000576) * t + 0.00200340) * t - 0.0001831) * t -
         46.836769) * t + 84381.406
    return ε
end

"""
    frame_tie(pos;direction)

Transforms a position vector from the dynamical reference system to ICRS, or vice versa.

# Argumenets
- `pos::AbstractVector`: Position vector, equitorial rectangular coordinates
- `direction::Symbol`: Transformation direction, `:dynamic2icrs` or `:icrs2dynamic`
"""
function frame_tie(pos::AbstractVector, direction::Symbol)
    @assert direction ∈ Set([:dynamic2icrs, :icrs2dynamic])
    if direction == :dynamic2icrs
        return frame_tie_rot * pos
    elseif direction == :icrs2dynamic
        return frame_tie_rot' * pos
    end
end

"""
    nutation(jd_tdb,pos)

Nutates equitorial rectangular coordinates from mean equator and equinox of epoch to true equator or
does the transformation in reverse given the flag `direction`.

# Arguments
- `jd_tdb::Real`: TDB Julian Date
- `pos::AbstractVector`: Position vector, geocentric equatorial rectangular coordinates, referred to mean equator and equinox of epoch

# Optional Arguments
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
- `direction::Symbol=:mean2true`: Perform the transformation in either the `:mean2true` or `:true2mean` direction
"""
function nutation(jd_tdb::Real, pos::AbstractVector; accuracy::Symbol=:full,
                  direction::Symbol=:mean2true)
    @assert direction ∈ Set([:mean2true, :true2mean])
    # Call e_tilt to get the obliquity and nutation angles
    oblm, oblt, _, psi, _ = e_tilt(jd_tdb; accuracy=accuracy)
    # Construct rotation matrix
    sobm, cobm = sincosd(oblm)
    sobt, cobt = sincosd(oblt)
    spsi, cpsi = sincos(psi * ASEC2RAD)
    rot = [cpsi -spsi*cobm -spsi*sobm
           spsi*cobt cpsi * cobm * cobt+sobm * sobt cpsi * sobm * cobt-cobm * sobt
           spsi*sobt cpsi * cobm * sobt-sobm * cobt cpsi * sobm * sobt+cobm * cobt]
    # Transform Position
    if direction == :mean2true
        return rot * pos
    elseif direction == :true2mean
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
@memoize function precession(jd_tdb1::Real, pos::AbstractVector, jd_tdb2::Real)
    @assert (jd_tdb1 == T0) || (jd_tdb2 == T0) "One of jd_tdb1, jd_tdb2 must be T0"
    eps0 = 84381.406
    # t is time in TDB centuries between the two epochs
    t = (jd_tdb2 - jd_tdb1) / 36525.0
    if jd_tdb2 == T0
        t *= -1
    end

    psia = ((((-0.0000000951 * t + 0.000132851) * t - 0.00114045) * t - 1.0790069) * t +
            5038.481507) * t
    omegaa = ((((+0.0000003337 * t - 0.000000467) * t - 0.00772503) * t + 0.0512623) * t -
              0.025754) * t + eps0
    chia = ((((-0.0000000560 * t + 0.000170663) * t - 0.00121197) * t - 2.3814292) * t +
            10.556403) * t

    eps0 = eps0 * ASEC2RAD
    psia = psia * ASEC2RAD
    omegaa = omegaa * ASEC2RAD
    chia = chia * ASEC2RAD

    sa, ca = sincos(eps0)
    sb, cb = sincos(-psia)
    sc, cc = sincos(-omegaa)
    sd, cd = sincos(chia)

    # Construct precession rotation matrix
    rot = [cd * cb-sb * sd * cc cd * sb * ca + sd * cc * cb * ca-sa * sd * sc cd*sb*sa+sd*cc*cb*sa+ca*sd*sc
           -sd * cb-sb * cd * cc -sd * sb * ca + cd * cc * cb * ca-sa * cd * sc -sd*sb*sa+cd*cc*cb*sa+ca*cd*sc
           sb*sc -sc * cb * ca-sa * cc -sc * cb * sa+cc * ca]

    # Perform rotation from epoch to J2000.0
    if jd_tdb2 == T0
        return rot' * pos
    else
        return rot * pos
    end
end

"""
    cel_pole(tjd)

Computes the celestial pole offsets for high-precision applications. Each set of
offsets is a correction to the modeled position of the pole for a specific date, derived
from observations and published by the IERS.

This function differs from the C version, where this returns the corrections instead of mutating global state.

# Arguments
- `tjd::Real`: TDB or TT Julian date for pole offsets
- `type::Symbol=:angular`: Type of pole offset. `:angular` for corrections to angular 
    coordinates of modeled pole referred to mean ecliptic of date, that is, 
    delta-delta-psi and delta-delta-epsilon. `:positional` for corrections to components
    of modeled pole unit vector referred to GCRS axes, that is, dx and dy.
- `dople1`: Value of celestial pole offset in first coordinate.
- `dpole2`: Value of celestial pole offset is second coordinate.
"""
function cel_pole(tjd::Real, type::Symbol, dpole1::Real, dpole2::Real)
    @assert type ∈ Set([:angular, :positional])
    if type == :angular
        psi_cor = dpole1 * 1e-3
        eps_cor = dpole2 * 1e-3
    elseif type == :positional
        # Components of modeled pole unit vector referred to GCRS axes
        dx = dpole1
        dy = dpole2
        t = (tjd - T0) / 36525.0
        # Compute sine_e of mean obliquity of date
        mean_ob = mean_obliq(tjd)
        sin_e = sin(mean_ob * ASEC2RAD)
        # Trivial model of pole trajectory is GCRS allows computation of dz
        x = (2004.190 * t) * ASEC2RAD
        dz = -(x + 0.5 * x * x * x) * dx
        # Form pole offset vector in GCRS
        dp1 = [dx, dy, dz] .* 1e-3 .* ASEC2RAD
        # Precess pole offset vector to mean equator and equinox of date
        dp2 = frame_tie(dp1, :icrs2dynamic)
        dp3 = precession(T0, dp2, tjd)
        # Compute delta-delta-psi and delta-delta-epsilon in arcseconds
        psi_cor = (dp3[1] / sin_e) / ASEC2RAD
        eps_cor = dp3[2] / ASEC3RAD
    end
    return psi_cor, eps_cor
end

"""
    e_tilt(jd_tdb)

Computes quantities related to the orientation of the Earth's rotation axis at Julian date `jd_tdb`

# Arguments
- `jd_tdb::Real`: TDB Julian Date

# Optional Arguments
- `accuracy::Symbol=:full`: Sets the accuracy level of `:full` or `:reduced`

# Returns
`(mobl,tobl,ee,dpsi,deps)` where
- `mobl`: Mean obliquity of the ecliptic in degrees
- `tobl`: True obliquity of the ecliptic in degrees
- `ee`: Equation of the equinoxes in seconds of time
- `dpsi`: Nutation in longitude on arcseconds
- `deps`: Nutation in obliquity in arcseconds
"""
@memoize function e_tilt(jd_tdb::Real; accuracy::Symbol=:full)
    # FIXME how do we include the PSI_COR and EPS_COR

    # Compute time in Julian centuries from epoch J2000.0
    t = (jd_tdb - T0) / 36525.0
    # Compute the nutation angles
    dp, de = nutation_angles(t; accuracy=accuracy)
    c_terms = ee_ct(jd_tdb, zero(typeof(jd_tdb)); accuracy=accuracy) / ASEC2RAD
    # Apply observed celestial pole offsets FIXME what do we do here?
    d_psi = dp # + PSI_COR
    d_eps = de # + EPS_COR
    # Compute mean obliquity of the ecliptic in arcseconds
    mean_ob = mean_obliq(jd_tdb)
    # Compute mean and true obliquity
    true_ob = mean_ob + d_eps
    # Convert obliquity values to degrees
    mean_ob /= 3600.0
    true_ob /= 3600.0
    # Compute equation of the equinoxes in seconds
    eq_eq = d_psi * cosd(mean_ob) + c_terms
    eq_eq /= 15.0
    # Return computed result
    return (mean_ob, true_ob, eq_eq, d_psi, d_eps)
end

"""
    ira_equinox(jd_tdb)

Computes the intermediate right ascension of the equinox at the input Julian date using an
analytical expression for the accumulated precession in right ascension. For the true equinox,
the result is the equation of the origins.

# Arguments
- `jd_tdb::Real`: TDB Julian Date

# Optional Arguments
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
- `equinox::Symbol=:mean`: Either `:mean` or `:true`
"""
@memoize function ira_equinox(jd_tdb::Real; accuracy::Symbol=:full,
                              equinox::Union{Symbol,Bool}=:mean)
    @assert equinox ∈ Set([:mean, :true])
    # Compute time in Julian centuries
    t = (jd_tdb - T0) / 36525.0
    # For the true equinox, obtain the equation of the equinoxes in seconds,
    # which includes the 'complementary terms'
    if equinox == :true
        _, _, eq_eq, _, _ = e_tilt(jd_tdb; accuracy=accuracy)
    elseif equinox == :mean
        eq_eq = 0.0
    end
    # Compute precession in RA in arcseconds taken from the reference
    prec_ra = 0.014506 +
              ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t +
               4612.156534) * t
    ra_eq = -(prec_ra / 15.0 + eq_eq) / 3600.0
    return ra_eq
end

"""
   cio_location(jd_tdb)

Computes the location of the celestial intermediate origin (CIO) for a given Julian date, as a
right ascension with respect to the true equinox of date.

# Arguments
- `jd_tdb::Real`: TDB Julian Date

# Optional Arguments
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
"""
@memoize function cio_location(jd_tdb::Real; accuracy::Symbol=:full)
    # NOTE: We are going to calculate the values instead of using the external dep for the time being
    # This results in ref_sys always being 2, true equator and equinox instead of interpolated GCRS
    ra_cio = -ira_equinox(jd_tdb; accuracy=accuracy, equinox=:true)
    return ra_cio
end

"""
    cio_basis(jd_tdb,ra_cio)

Computes the orthonomal basis vectors, with respect to the GCRS of the celestial intermediate system defined
by the celestial intermediate pole and origin.

# Arguments
- `jd_tdb::Real`: TDB Julian Date
- `ra_cio::Real`: Right ascension of the CIO at epoch (hours)

# Optional Arguments
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
"""
@memoize function cio_basis(jd_tdb::Real, ra_cio::Real; accuracy::Symbol=:full)
    # Compute unit vector z towards celestial pole
    z0 = [0.0, 0.0, 1.0]
    w1 = nutation(jd_tdb, z0; accuracy=accuracy, direction=:true2mean)
    w2 = precession(jd_tdb, w1, T0)
    z = frame_tie(w2, :dynamic2icrs)
    # Compute unit vectors x and y
    # First construct unit vector towards CIO in equator and equinox of date system
    w0 = [cosd(ra_cio * 15), sind(ra_cio * 15), 0]
    # Rotate the vector into the GCRS to form unit vector x
    w1 = nutation(jd_tdb, w0; accuracy=accuracy, direction=:true2mean)
    w2 = precession(jd_tdb, w1, T0)
    x = frame_tie(w2, :dynamic2icrs)
    # Compute unit vector y orthogonal to x and z
    y = z × x
    return x, y, z
end

"""
    era(jd_high)

This function returns the value of the earth rotation angle (theta in degrees) for a given UT1 Julian date.

# Arguments
- `jd_high::Real`: High-order part of UT1 Julian date
- `jd_low::Real=0.0`: Low-order part of UT1 Julian date
"""
function era(jd_high::Real, jd_low::Real=0.0)
    thet1 = 0.7790572732640 + 0.00273781191135448 * (jd_high - T0)
    thet2 = 0.00273781191135448 * jd_low
    thet3 = rem(jd_high, 1.0) + rem(jd_low, 1.0)
    theta = rem(thet1 + thet2 + thet3, 1.0) * 360.0
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

Computes the Greenwich sidereal time, either mean or apparent, at Julian date `jd_high` + `jd_low` in hours.
See Chapter 5 in the NOVASC manual for information between equinox and CIO-based methods.

# Arguments
- `jd_high::Real`: high-order part of UT1 Julian date
- `jd_low::Real=0`: low-order part of UT1 Julian date
- `delta_t::Real=0`: Difference TT-UT1 at `jd_high` + `jd_low` in seconds

# Optional arguments
- `gst_type::Symbol=:mean`: Return results as mean (`:mean`) or apparent (`:apparent`) time
- `method::Symbol=:CIO`: Computation method, `:CIO`-based or `:equinox`-based
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
"""
@memoize function sidereal_time(jd_high::Real, jd_low::Real=0.0, delta_t::Real=0.0;
                                gst_type::Symbol=:mean, method::Symbol=:CIO,
                                accuracy::Symbol=:full)
    @assert method ∈ Set([:equinox, :CIO])
    @assert gst_type ∈ Set([:mean, :apparent])

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
    if ((gst_type == :mean) && (method == :CIO)) ||
       ((gst_type == :apparent) && (method == :equinox))
        _, _, ee, _, _ = e_tilt(jd_tdb; accuracy=accuracy)
        eqeq = ee * 15
    else
        eqeq = 0
    end

    if method == :CIO
        # Obtain basis vectors, in the GCRS of the celestial intermediate system
        ra_cio = cio_location(jd_tdb; accuracy=accuracy)
        x, y, z = cio_basis(jd_tdb, ra_cio; accuracy=accuracy)
        # Compute the direction of the true equinox in the GCRS
        w1 = nutation(jd_tdb, [1, 0, 0]; direction=:true2mean, accuracy=accuracy)
        w2 = precession(jd_tdb, w1, T0)
        eq = frame_tie(w2, :dynamic2icrs)
        # Compute the hour angle of the equinox w.r.t. the TIO meridian
        ha_eq = theta - atand(eq ⋅ y, eq ⋅ x)
        # For mean, subtract the equation of the equinoxes
        ha_eq -= eqeq / 240
        ha_eq = rem(ha_eq, 360) / 15
        if ha_eq < 0
            ha_eq += 24
        end
        gst = ha_eq
    elseif method == :equinox
        # Precession-in-RA terms in mean sidereal time with coefficients in arcseconds
        st = eqeq +
             0.014506 +
             ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t +
              4612.156534) * t
        gst = rem((st / 3600.0 + theta), 360.0) / 15.0
        if gst < 0.0
            gst += 24.0
        end
    end
    return gst
end

"""
    wobble(tjd,xp,yp,pos)

Corrects a vector in the ITRS (rotating Earth-fixed system)
for polar motion, and also corrects the longitude origin
(by a tiny amount) to the Terrestrial Intermediate Origin
(TIO).  The ITRS vector is thereby transformed to the terrestrial
intermediate system, based on the true (rotational) equator and
TIO.  Because the true equator is the plane orthogonal to the
direction of the Celestial Intermediate Pole (CIP), the components
of the output vector are referred to z and x axes toward the CIP
and TIO, respectively.

# Arguments
- `tjd::Real`: TT or UT1 Julian Date
- `xp::Real`: Conventionally-defined X coordinate of the CIP in arcseconds
- `yp::Real`: Conventionally-defined Y coordinate of the CIP in arcseconds
- `pos::AbstractVector`: Position vector, geocentric equatorial rectangular coordinates

# Optional Arguments
- `direction::Symbol=:itrs2terr`: Either `:itrs2terr` or `:terr2itrs`
"""
function wobble(tjd::Real, xp::Real, yp::Real, pos::AbstractVector;
                direction::Symbol=:itrs2terr)
    @assert direction ∈ Set([:itrs2terr, :terr2itrs])
    xpole = xp * ASEC2RAD
    ypole = yp * ASEC2RAD
    # Compute approximate longitude of the TIO
    t = (tjd - T0) / 36525.0
    sprime = -47.0e-6 * t
    tiolon = -sprime * ASEC2RAD
    # Compute elements of rotation matrix
    sinx, cosx = sincos(xpole)
    siny, cosy = sincos(ypole)
    sinl, cosl = sincos(tiolon)

    xx = cosx * cosl
    yx = sinx * siny * cosl + cosy * sinl
    zx = -sinx * cosy * cosl + siny * sinl
    xy = -cosx * sinl
    yy = -sinx * siny * sinl + cosy * cosl
    zy = sinx * cosy * sinl + siny * cosl
    xz = sinx
    yz = -cosx * siny
    zz = cosx * cosy

    rot_mat = @SMatrix [xx yx zx; xy yy zy; xz yz zz]

    # Perform rotation
    if direction == :itrs2terr
        return rot_mat * pos
    elseif direction == :terr2itrs
        return rot_mat' * pos
    end
end

"""
    spin(angle,pos)

Transfroms a vector from one coordiate system to another with the same origin and axes rotated about the z-axis.

# Arguments
- `angle::Real`: Angle of rotation, positive counterclockwise when viewed from +z, in degrees
- `pos::AbstractVector`: Position Vector
"""
function spin(angle::Real, pos::AbstractVector)
    sinang, cosang = sincosd(angle)
    xx = cosang
    yx = sinang
    zx = 0.0
    xy = -sinang
    yy = cosang
    zy = 0.0
    xz = 0.0
    yz = 0.0
    zz = 1.0
    rot_mat = @SMatrix [xx yx zx; xy yy zy; xz yz zz]
    return rot_mat * pos
end

"""
    ter2cel(jd_ut_high,jd_ut_low,delta_t,vec)

Rotates a vector from the terrestrial to the celestial system. Specifically, it transforms
a vector in ITRS to GCRS by applying rotations for polar motion, Earth rotation, nutation, precession, and
the dynamical-to-GCRS frame tie.

# Arguments
- `jd_ut_high::Real`: High order part of UT1 Julian date
- `jd_ut_low::Real`: Low order part of UT1 Julian date
- `delta_t::Real`: Value of TT-UT1 at the input date
- `vec::AbstractVector`: Position vector, geocentric equatorial rectangular coordinates, reference axes set by `option`

# Optional Arguments
- `method::Symbol=:CIO`: Computation method, `:CIO`-based or `:equinox`-based
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
- `option::Symbol=:GRCS`: Sets position vector to be referred to `:GCRS` axes or `:equinox` of date and equator
- `xp::Real=0`: x coordinate of the celestial intermediate pole in arcseconds
- `yp::Real=0`: y coordinate of the celestial intermediate pole in arcseconds
"""
function ter2cel(jd_ut_high::Real, jd_ut_low::Real, delta_t::Real, vec::AbstractVector;
                 method::Symbol=:CIO, accuracy::Symbol=:full, option::Symbol=:GCRS,
                 xp::Real=0, yp::Real=0)
    @assert method ∈ Set([:CIO, :equinox])
    @assert option ∈ Set([:GCRS, :equinox])
    # Compute the TT Julian date
    jd_ut1 = jd_ut_high + jd_ut_low
    jd_tt = jd_ut1 + (delta_t / 86400.0)
    # Compute the TDB Julian date corresponding to the input UT1 Julian date
    jd_tdb = jd_tt
    _, secdiff = tdb2tt(jd_tdb)
    jd_tdb = jd_tt + secdiff / 86400.0
    # Apply polar motion
    if iszero(xp) && iszero(yp)
        v1 = vec
    else
        v1 = wobble(jd_tdb, xp, yp, vec; direction=:itrs2terr)
    end
    if method == :CIO
        # Obtain basis vectors, in the GCRS of the celestial intermediate system
        r_cio = cio_location(jd_tdb; accuracy=accuracy)
        x, y, z = cio_basis(jd_tdb, r_cio; accuracy=accuracy)
        # Compute and apply earth rotation angle, `theta` transforming the vector to the celestial intermediate system
        θ = era(jd_ut_high, jd_ut_low)
        v2 = spin(-θ, v1)
        # Transform the vector from the celestial intermediate system to the GCRS
        vec2 = hcat(x, y, z) * v2
    elseif method == :equinox
        gast = sidereal_time(jd_ut_high, jd_ut_low, delta_t; gst_type=:apparent,
                             method=:equinox, accuracy=accuracy)
        v2 = spin(-gast * 15, v1)
        if option == :equinox
            vec2 = v2
        elseif option == :GCRS
            # Apply precession, nutation, frame tie
            v3 = nutation(jd_tdb, v2; accuracy=accuracy, direction=:true2mean)
            v4 = precession(jd_tdb, v3, T0)
            vec2 = frame_tie(v4, :dynamic2icrs)
        end
    end
    return vec2
end

"""
    OneSurface(lat,lon,h,temp,p)

Structure holding the observer's location on the surface of the Earth. 
Atmospheric parameters are optional and are only used by `equ2hor` in the
calculation of refraction.

# Arguments
- `latitude::Real`: ITRS latitude in degrees; north positive
- `longitude::Real`: ITRS longitude in degrees; east positive
- `height::Real`: Height of observer (meters)
- `temperature::Real=0`: Temperature (degrees Celsius)
- `pressure::Real=0`: Atmospheric pressure (millibars)
"""
struct OnSurface{T<:Real}
    latitude::T
    longitude::T
    height::T
    temperature::T
    pressure::T
end

OnSurface(lat::T, lon::T, h::T) where {T<:Real} = OnSurface(lat, lon, h, zero(T), zero(T))

"""
    refract(location,zd_obs)

Computes the atmospheric refraction in zenith distance.
This is approximate for optical wavelengths.

# Arguments
- `location::OnSurface`: The observation location
- `zd_obs::Real`: Observed zenith distance, in degrees

# Optional Arguments
- `ref_option::Symbol=:standard`: Either `:standard` to use standard atmospheric conditions, or `:location` to use the conditions included in `location`.
"""
function refract(location::OnSurface, zd_obs::Real; ref_option::Symbol=:standard)
    @assert ref_option ∈ Set([:standard, :location])
    # We only care about zenith distances between 0.1 and 91 degrees
    if 0.1 <= zd_obs <= 91
        if ref_option == :location
            p = location.pressure
            t = location.temperature
        elseif ref_option == :standard
            p = 1010.0 * exp(-location.height / 9.1e3)
            t = 10.0
        end
        h = 90.0 - zd_obs
        r = 0.016667 / tand(h + 7.31 / (h + 4.4))
        refr = r * (0.28 * p / (t + 273.0))
    else
        refr = 0.0
    end
    return refr
end

"""
    equ2hor(jd_ut1,delta_t,ra,dec)

This function transforms topocentric right ascension and
declination to zenith distance and azimuth.  It uses a method
that properly accounts for polar motion, which is significant at
the sub-arcsecond level.  This function can also adjust
coordinates for atmospheric refraction.

`xp` and `yp` are only needed if sub-arcsecond accuracy is needed.

# Arguments
- `jd_ut1::Real`: UT1 Julian Date
- `delta_t::Real`: Difference TT-UT1 ast `jd_ut1` in seconds
- `ra::Real`: Topocentric right ascension in hours, referred to true equator
- `dec::Real`: Topocentric declination in degrees, referred to true equatior
- `location::OnSurface`: Structure containing observer's location details

# Optional arguments
- `accuracy::Symbol=:full`: Either `:full` or `:reduced` accuracy
- `xp::Real=0`: x coordinate of the celestial intermediate pole in arcseconds
- `yp::Real=0`: y coordinate of the celestial intermediate pole in arcseconds
- `ref_option::Symbol=:none`: Whether to include refraction in the calculations.
Either `:none`, `:standard` to use standard atmospheric conditions, or `:location` to use the conditions
included in `location`.

# Returns
`(zd,az,rar,decr)` where
- `zd`: Topocentric zenith distance in degrees
- `az`: Topocentric azimuth (east of north) in degrees
- `rar`: Topocentric RA of object of interest in hours
- `decr`: Topocentric declination of object of interest in degrees
"""
function equ2hor(jd_ut1::Real, delta_t::Real, ra::Real, dec::Real, location::OnSurface;
                 accuracy::Symbol=:full, xp::Real=0.0, yp::Real=0.0,
                 ref_option::Symbol=:none)
    @assert ref_option ∈ Set([:none, :standard, :location])
    # Initialize things
    refr = 0.0
    rar = ra
    decr = dec
    # Trig it up
    sinlat, coslat = sincosd(location.latitude)
    sinlon, coslon = sincosd(location.longitude)
    sindc, cosdc = sincosd(dec)
    sinra, cosra = sincosd(ra * 15)
    # Set up orthonomarl basis vectors in local Earth-fixed system
    uze = [coslat * coslon, coslat * sinlon, sinlat]
    # Define vector towards local north in Earth-fixed system
    une = [-sinlat * coslon, -sinlat * sinlon, coslat]
    # Define vector towards local west in Earth-fixed system
    uwe = [sinlon, -coslon, 0.0]
    # Obtain vectors in celestial system and rotate earth-fixed orthonormal basis to celestial system
    uz = ter2cel(jd_ut1, 0.0, delta_t, uze; accuracy=accuracy, xp=xp, yp=yp,
                 option=:equinox, method=:equinox)
    un = ter2cel(jd_ut1, 0.0, delta_t, une; accuracy=accuracy, xp=xp, yp=yp,
                 option=:equinox, method=:equinox)
    uw = ter2cel(jd_ut1, 0.0, delta_t, uwe; accuracy=accuracy, xp=xp, yp=yp,
                 option=:equinox, method=:equinox)
    # Define unit vector `p` towards object in celestial system w.r.t equator and equinox of date
    p = [cosdc * cosra, cosdc * sinra, sindc]
    # Compute coordinates of object w.r.t orthonomral basis
    # Compuute projected components of `p` onto rotated Earth-fixed basis vectors
    pz = p[1] * uz[1] + p[2] * uz[2] + p[3] * uz[3]
    pn = p[1] * un[1] + p[2] * un[2] + p[3] * un[3]
    pw = p[1] * uw[1] + p[2] * uw[2] + p[3] * uw[3]
    # Compute azimuth and zenith distance
    proj = sqrt(pn^2 + pw^2)
    if proj > 0
        az = -atand(pw, pn)
    end
    if az < 0
        az += 360
    end
    if az >= 360
        az -= 360
    end
    zd = atand(proj, pz)
    # Apply atmospheric refraction if requested
    if ref_option != :none
        # Get refraction in zenith distance
        # This is iterative because refraction algorithms are always a function of observed zenith
        zd0 = zd
        while true
            zd1 = zd
            refr = refract(location, zd; ref_option=ref_option)
            zd = zd0 - refr
            abs(zd - zd1) > 3.0e-5 || break
        end
        # Apply refraction to celestial coordinates of object
        if (refr > 0) && (zd > 3.0e-4)
            # Shift position vector of object in celestial system to account for refraction
            sinzd, coszd = sincosd(zd)
            sinzd0, coszd0 = sincosd(zd0)
            # Compute refracted position vector
            pr = @. ((p - coszd0 * uz) / sinzd0) * sinzd + uz * coszd
            # Compute refracted ra and dec
            proj = sqrt(pr[1]^2 + pr[2]^2)
            if proj > 0.0
                rar = atand(pr[2], pr[1]) / 15.0
            end
            if rar < 0.0
                rar += 24.0
            end
            if rar >= 24.0
                rar -= 24.0
            end
            decr = atand(pr[3], proj)
        end
    end
    return zd, az, rar, decr
end
