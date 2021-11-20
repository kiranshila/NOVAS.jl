var documenterSearchIndex = {"docs":
[{"location":"api/#API-(Exported-Functions)","page":"API","title":"API (Exported Functions)","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [NOVAS]\nOrder = [:function,:type]\nPrivate = false","category":"page"},{"location":"api/#NOVAS.cio_basis-Tuple{Real, Real}","page":"API","title":"NOVAS.cio_basis","text":"cio_basis(jd_tdb,ra_cio)\n\nComputes the orthonomal basis vectors, with respect to the GCRS of the celestial intermediate system defined by the celestial intermediate pole and origin.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\nra_cio::Real: Right ascension of the CIO at epoch (hours)\n\nOptional Arguments\n\naccuracy::Symbol=:full: Either :full or :reduced accuracy\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.cio_location-Tuple{Real}","page":"API","title":"NOVAS.cio_location","text":"ciolocation(jdtdb)\n\nComputes the location of the celestial intermediate origin (CIO) for a given Julian date, as a right ascension with respect to the true equinox of date.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\n\nOptional Arguments\n\naccuracy::Symbol=:full: Either :full or :reduced accuracy\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.e_tilt-Tuple{Real}","page":"API","title":"NOVAS.e_tilt","text":"e_tilt(jd_tdb)\n\nComputes quantities related to the orientation of the Earth's rotation axis at Julian date jd_tdb\n\nArguments\n\njd_tdb::Real: TDB Julian Date\n\nOptional Arguments\n\naccuracy::Symbol=:full: Sets the accuracy level of :full or :reduced\n\nReturns\n\n(mobl,tobl,ee,dpsi,deps) where\n\nmobl: Mean obliquity of the ecliptic in degrees\ntobl: True obliquity of the ecliptic in degrees\nee: Equation of the equinoxes in seconds of time\ndpsi: Nutation in longitude on arcseconds\ndeps: Nutation in obliquity in arcseconds\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.ee_ct-Union{Tuple{T}, Tuple{T, T}} where T<:Real","page":"API","title":"NOVAS.ee_ct","text":"ee_ct(jd)\n\nComputes the \"complementary terms\" of the equation of the equinoxes.\n\nArguments\n\njd_high::Real: High-order part of UT1 Julian date\njd_low::Real=0.0: Low-order part of UT1 Julian date\n\nOptional Arguments\n\naccuracy::Symbol=:full: Sets the accuracy level of :full or :reduced\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.era","page":"API","title":"NOVAS.era","text":"era(jd_high)\n\nThis function returns the value of the earth rotation angle (theta in degrees) for a given UT1 Julian date.\n\nArguments\n\njd_high::Real: High-order part of UT1 Julian date\njd_low::Real=0.0: Low-order part of UT1 Julian date\n\n\n\n\n\n","category":"function"},{"location":"api/#NOVAS.frame_tie-Tuple{AbstractVector{T} where T, Symbol}","page":"API","title":"NOVAS.frame_tie","text":"frame_tie(pos;direction)\n\nTransforms a position vector from the dynamical reference system to ICRS, or vice versa.\n\nArgumenets\n\npos::AbstractVector: Position vector, equitorial rectangular coordinates\ndirection::Symbol: Transformation direction, :dynamic2icrs or :icrs2dynamic\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.fund_args-Tuple{Real}","page":"API","title":"NOVAS.fund_args","text":"fund_args(t)\n\nComputes the fundamental arguments (mean elements) of the Sun and Moon.\n\nArguments\n\n-t::Real: TDB time in Julian centuries since J2000.0\n\nReturns\n\n[l,l′,F,D,Ω] where -l: Mean anomaly of the Moon -l′: Mean anomaly of the Sun -F: mean argument of the latitude of the Moon -D: Mean elongation of the Moon from the Sun -Ω Mean longitude of the Moon's ascending node\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.iau2000a","page":"API","title":"NOVAS.iau2000a","text":"iau2000a(jd)\n\nCompute the forced nutation of the non-rigid earth based on the IAU 2000A nutation model.\n\nArguments\n\n-jd_high::Real: High order part of the TT Julian date -jd_low::Real=0.0: Low order part of the TT Julian date\n\nReturn\n\n(dpsi,deps) where -dpsi: Nutation (luni-solar + planetary) in longitude, in radians -deps: Nutation (luni-solar + planetary) in obliquity, in radians\n\n\n\n\n\n","category":"function"},{"location":"api/#NOVAS.ira_equinox-Tuple{Real}","page":"API","title":"NOVAS.ira_equinox","text":"ira_equinox(jd_tdb)\n\nComputes the intermediate right ascension of the equinox at the input Julian date using an analytical expression for the accumulated precession in right ascension. For the true equinox, the result is the equation of the origins.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\n\nOptional Arguments\n\naccuracy::Symbol=:full: Either :full or :reduced accuracy\nequinox::Symbol=:mean: Either :mean or :true\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.mean_obliq-Tuple{Real}","page":"API","title":"NOVAS.mean_obliq","text":"mean_obliq(jd_tdb)\n\nComputes the mean obliquity of the ecliptic\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.nu2000k","page":"API","title":"NOVAS.nu2000k","text":"nu2000k(jd)\n\nCompute the forced nutation of the non-rigid earth based on the NU2000K nutation model.\n\nArguments\n\n-jd_high::Real: High order part of the TT Julian date -jd_low::Real=0.0: Low order part of the TT Julian date\n\nReturn\n\n(dpsi,deps) where -dpsi: Nutation (luni-solar + planetary) in longitude, in radians -deps: Nutation (luni-solar + planetary) in obliquity, in radians\n\n\n\n\n\n","category":"function"},{"location":"api/#NOVAS.nutation-Tuple{Real, AbstractVector{T} where T}","page":"API","title":"NOVAS.nutation","text":"nutation(jd_tdb,pos)\n\nNutates equitorial rectangular coordinates from mean equator and equinox of epoch to true equator or does the transformation in reverse given the flag direction.\n\nArguments\n\njd_tdb::Real: TDB Julian Date\npos::AbstractVector: Position vector, geocentric equatorial rectangular coordinates, referred to mean equator and equinox of epoch\n\nOptional Arguments\n\naccuracy::Symbol=:full: Either :full or :reduced accuracy\ndirection::Symbol=:mean2true: Perform the transformation in either the :mean2true or :true2mean direction\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.nutation_angles-Tuple{Real}","page":"API","title":"NOVAS.nutation_angles","text":"nutation_angles(t)\n\nComputes the values for nutation in longitude and nutation in obliquity for a given TDB julian date.\n\nAruments\n\nt::Real: TDB time in Julian centuries since J2000.0\n\nOptional Arguments\n\naccuracy::Symbol=:full: Sets the accuracy level of :full or :reduced\n\nReturns\n\n(dpsi,deps) where\n\ndpsi: Nutation in longitude in arcseconds\ndeps: Nutation in obliquity in arcseconds\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.precession-Tuple{Real, AbstractVector{T} where T, Real}","page":"API","title":"NOVAS.precession","text":"precession(jd_tdb1,pos,jd_tdb2)\n\nPrecesses equatorial rectangular coordinates from one epoch to another. One of the two epochs must be J200.0.\n\nArguments\n\njd_tdb1::Real: TDB Julian Date of first epoch\njd_tdb2::Real: TDB Julian Date of second epoch\npos::AbstractVector: Position vector, geocentric equatorial rectangular coordinates, referred to mean dynamical equator and equinox of first epoch\n\n\n\n\n\n","category":"method"},{"location":"api/#NOVAS.sidereal_time","page":"API","title":"NOVAS.sidereal_time","text":"sidereal_time(jd_ut1,0,delta_t)\n\nComputes the Greenwich sidereal time, either mean or apparent, at Julian date jd_high + jd_low in hours\n\nArguments\n\njd_high::Real: high-order part of UT1 Julian date\njd_low::Real=0: low-order part of UT1 Julian date\ndelta_t::Real=0: Difference TT-UT1 at jd_high + jd_low in seconds\n\nOptional arguments\n\ngst_type::Symbol=:mean: Return results as mean (:mean) or apparent (:apparent) time\nmethod::Symbol=:CIO: Computation method, CIO-based (:CIO) or equinox-based (:equinox)\naccuracy::Symbol=:full: Either :full or :reduced accuracy\n\n\n\n\n\n","category":"function"},{"location":"api/#NOVAS.tdb2tt-Tuple{Real}","page":"API","title":"NOVAS.tdb2tt","text":"tdb2tt(tdb_jd)\n\nComputes the Terrestrial Time (TT) or Terrestrial Dynamic Time (TDT) Julian date corresponding to a Barrycentric Dynamical Time (TDB) Julian Date)\n\nArguments\n\ntdb_jd::Real: TDB Julian date\n\n\n\n\n\n","category":"method"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Calculating-Sidereal-Time","page":"Examples","title":"Calculating Sidereal Time","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example, we will calculate the Greenwich mean sidereal time.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using NOVAS","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"First, we need to pull in AstroTime.jl to give us the most preceise time correction data and transformations. We will also call it's update routine.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using AstroTime\nAstroTime.update()","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Next, we will use the current Julian UT1 date.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"ut1_now = from_utc(now(); scale=UT1)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Then, we will use AstroTime to convert this date to a Julian date and give us the current offset to Terrestrial Time","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"jd_high, jd_low = julian_twopart(ut1_now) .|> value\nΔT = getoffset(ut1_now,TT)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Finally, we will call sidereal_time to find the Greenwich mean sidereal time using the CIO method with full accuracy","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"sidereal_time(jd_high,jd_low,ΔT)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = NOVAS","category":"page"},{"location":"#NOVAS","page":"Home","title":"NOVAS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NOVAS.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"examples.md\", \"api.md\"]\nDepth = 3","category":"page"}]
}