var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = NOVAS","category":"page"},{"location":"#NOVAS","page":"Home","title":"NOVAS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NOVAS.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NOVAS]\nOrder = [:function,:type]\nPrivate = false","category":"page"},{"location":"#NOVAS.cio_location-Tuple{Real}","page":"Home","title":"NOVAS.cio_location","text":"ciolocation(jdtdb)\n\nComputes the location of the celestial intermediate origin (CIO) for a given Julian date, as a right ascension with respect to the true equinox of date.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\n\nOptional Arguments\n\naccuracy::Accuracy=full: Either full or reduced accuracy\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.e_tilt-Tuple{Real}","page":"Home","title":"NOVAS.e_tilt","text":"e_tilt(jd_tdb)\n\nComputes quantities related to the orientation of the Earth's rotation axis at Julian date jd_tdb\n\nArguments\n\njd_tdb::Real: TDB Julian Date\n\nOptional Arguments\n\naccuracy::Accuracy=full: Sets the accuracy level of full or reduced\n\nReturns\n\n(mobl,tobl,ee,dpsi,deps) where\n\nmobl: Mean obliquity of the ecliptic in degrees\ntobl: True obliquity of the ecliptic in degrees\nee: Equation of the equinoxes in seconds of time\ndpsi: Nutation in longitude on arcseconds\ndeps: Nutation in obliquity in arcseconds\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.era-Tuple{Real, Real}","page":"Home","title":"NOVAS.era","text":"era(jd_high)\n\nThis function returns the value of the earth rotation angle (theta in degrees) for a given UT1 Julian date.\n\nArguments\n\njd_high::Real: High-order part of UT1 Julian date\njd_low::Real: Low-order part of UT1 Julian date\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.frame_tie","page":"Home","title":"NOVAS.frame_tie","text":"frame_tie(pos;direction)\n\nTransforms a position vector from the dynamical reference system to ICRS, or vice versa.\n\nArgumenets\n\npos::AbstractVector: Position vector, equitorial rectangular coordinates\ndirection::FrameTieDirection=dynamic2icrs\n\n\n\n\n\n","category":"function"},{"location":"#NOVAS.fund_args-Tuple{Real}","page":"Home","title":"NOVAS.fund_args","text":"fund_args(t)\n\nComputes the fundamental arguments (mean elements) of the Sun and Moon.\n\nArguments\n\n-t::Real: TDB time in Julian centuries since J2000.0\n\nReturns\n\n[l,l′,F,D,Ω] where -l: Mean anomaly of the Moon -l′: Mean anomaly of the Sun -F: mean argument of the latitude of the Moon -D: Mean elongation of the Moon from the Sun -Ω Mean longitude of the Moon's ascending node\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.iau2000a","page":"Home","title":"NOVAS.iau2000a","text":"iau2000a(jd)\n\nCompute the forced nutation of the non-rigid earth based on the IAU 2000A nutation model.\n\nArguments\n\n-jd_high::Real: High order part of the TT Julian date -jd_low::Real=0.0: Low order part of the TT Julian date\n\nReturn\n\n(dpsi,deps) where -dpsi: Nutation (luni-solar + planetary) in longitude, in radians -deps: Nutation (luni-solar + planetary) in obliquity, in radians\n\n\n\n\n\n","category":"function"},{"location":"#NOVAS.ira_equinox-Tuple{Real}","page":"Home","title":"NOVAS.ira_equinox","text":"ira_equinox(jd_tdb)\n\nComputes the intermediate right ascension of the equinox at the input Julian date using an analytical expression for the accumulated precession in right ascension. For the true equinox, the result is the equation of the origins.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\n\nOptional Arguments\n\naccuracy::Accuracy=full: Either full or reduced accuracy\nequinox::EquinoxType=mean_equinox: Either mean_equinox or true_equinox\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.nutation-Tuple{Real, AbstractVector{T} where T}","page":"Home","title":"NOVAS.nutation","text":"nutation(jd_tdb)\n\nNutates equitorial rectangular coordinates from mean equator and equinox of epoch to true equator or does the transformation in reverse given the flag direction.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\npos::AbstractVector: Position vector, geocentric equatorial rectangular coordinates, referred to mean equator and equinox of epoch\n\nOptional Arguments\n\naccuracy::Accuracy=full: Either full or reduced accuracy\ndirection::NutationDirection=mean2true: Perform the transformation in either the mean2true or true2mean direction\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.nutation_angles-Tuple{Real}","page":"Home","title":"NOVAS.nutation_angles","text":"nutation_angles(t)\n\nComputes the values for nutation in longitude and nutation in obliquity for a given TDB julian date.\n\nAruments\n\nt::Real: TDB time in Julian centuries since J2000.0\n\nOptional Arguments\n\naccuracy::Accuracy=full: Sets the accuracy level of full or reduced\n\nReturns\n\n(dpsi,deps) where\n\ndpsi: Nutation in longitude in arcseconds\ndeps: Nutation in obliquity in arcseconds\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.precession-Tuple{Real, AbstractVector{T} where T, Real}","page":"Home","title":"NOVAS.precession","text":"precession(jd_tdb1,pos,jd_tdb2)\n\nPrecesses equatorial rectangular coordinates from one epoch to another. One of the two epochs must be J200.0.\n\nArguments\n\njd_tdb1::Real: TDB Julian Date of first epoch\njd_tdb2::Real: TDB Julian Date of second epoch\npos::AbstractVector: Position vector, geocentric equatorial rectangular coordinates, referred to mean dynamical equator and equinox of first epoch\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.read_iau2000a-Tuple{}","page":"Home","title":"NOVAS.read_iau2000a","text":"read_iau2000a()\n\nReads the included iau2000a nutation model into a DataFrame and returns them as (planetary_nutations,lunisolar_nutations). This is meant to be a high-level interface for working with the iau2000a data.\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.sidereal_time-Tuple{Real, Real, Real}","page":"Home","title":"NOVAS.sidereal_time","text":"sidereal_time(jd_ut1,0,delta_t)\n\nComputes the Greenwich sidereal time, either mean or apparent, at Julian date jd_high + jd_low in hours\n\nArguments\n\njd_high::Real: high-order part of UT1 Julian date\njd_low::Real: low-order part of UT1 Julian date\ndelta_t::Real: Difference TT-UT1 at jd_high + jd_low in seconds\n\nOptional arguments\n\ngst_type::GSTType=mean_gst: Return results as mean (mean_gst) or apparent (apparent_gst) time\nmethod::SiderealMethod=equinox: Computation method, CIO-based (CIO) or equinox-based (equinox)\naccuracy::Accuracy=full: Either full or reduced accuracy\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.tdb2tt-Tuple{Real}","page":"Home","title":"NOVAS.tdb2tt","text":"tdb2tt(tdb_jd)\n\nComputes the Terrestrial Time (TT) or Terrestrial Dynamic Time (TDT) Julian date corresponding to a Barrycentric Dynamical Time (TDB) Julian Date)\n\nArguments\n\ntdb_jd::Real: TDB Julian date\n\n\n\n\n\n","category":"method"}]
}
