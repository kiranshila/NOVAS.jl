# Barycentric Dynamic Time Julian date of epoch J2000.0
const T0 = 2451545.0
# Speed of light in meters/second is a defining physical constant.
const C = 299792458.0
# Light-time for one astronomical unit (AU) in seconds, from DE-405.
const AU_SEC = 499.0047838061
# Speed of light in AU/day.  Value is 86400 / AU_SEC.
const C_AUDAY = 86400/AU_SEC
# Astronomical unit in meters.  Value is AU_SEC * C.
const AU = AU_SEC*C
# Astronomical Unit in kilometers.
const AU_KM = 1.4959787069098932e+8
# Heliocentric gravitational constant in meters^3 / second^2, from
# DE-405.
const GS = 1.32712440017987e+20
# Geocentric gravitational constant in meters^3 / second^2, from
# DE-405.
const GE = 3.98600433e+14
# Radius of Earth in meters from IERS Conventions (2003).
const ERAD = 6378136.6
# Earth ellipsoid flattening from IERS Conventions (2003).
# Value is 1 / 298.25642.
const F = 1/298.25642
# Rotational angular velocity of Earth in radians/sec from IERS
# Conventions (2003).
const ANGVEL = 7.2921150e-5
# Reciprocal masses of solar system bodies, from DE-405
# (Sun mass / body mass).
# MASS[0] = Earth/Moon barycenter, MASS[1] = Mercury, ...,
# MASS[9] = Pluto, MASS[10] = Sun, MASS[11] = Moon.
const RMASS = [328900.561400, 6023600.0, 408523.71,
               332946.050895, 3098708.0, 1047.3486, 3497.898, 22902.98,
               19412.24, 135200000.0, 1.0, 27068700.387534]
# Number of arcseconds in 360 degrees.
const ASEC360 = 1296000.0
# Angle conversion constants
const ASEC2RAD = 4.848136811095359935899141e-6
const DEG2RAD = 0.017453292519943296
const RAD2DEG = 57.295779513082321
