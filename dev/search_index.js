var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = NOVAS","category":"page"},{"location":"#NOVAS","page":"Home","title":"NOVAS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NOVAS.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NOVAS]\nOrder = [:function,:type]\nPrivate = false","category":"page"},{"location":"#NOVAS.e_tilt-Tuple{Real}","page":"Home","title":"NOVAS.e_tilt","text":"e_tilt!(jd_tdb)\n\nComputes quantities related to the orientation of the Earth's rotation axis at Julian date jd_tdb\n\nArguments\n\njd_tdb::Real: TDB Julian Date\naccuracy=:full: Sets the accuracy level of :full or :reduced\n\nReturns\n\n(mobl,tobl,ee,dpsi,deps) where\n\nmobl: Mean obliquity of the ecliptic in degrees\ntobl: True obliquity of the ecliptic in degrees\nee: Equation of the equinoxes in seconds of time\ndpsi: Nutation in longitude on arcseconds\ndeps: Nutation in obliquity in arcseconds\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.era-Tuple{Real, Real}","page":"Home","title":"NOVAS.era","text":"era(jd_high)\n\nThis function returns the value of the earth rotation angle (theta in degrees) for a given UT1 Julian date.\n\nArguments\n\njd_high::Real: High-order part of UT1 Julian date\njd_low::Real: Low-order part of UT1 Julian date\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.sidereal_time-Tuple{Real, Real, Real}","page":"Home","title":"NOVAS.sidereal_time","text":"sidereal_time(jd_ut1,0,delta_t)\n\nComputes the Greenwich sidereal time, either mean or apparent, at Julian date jd_high + jd_low in hours\n\nArguments\n\njd_high::Real: high-order part of UT1 Julian date\njd_low::Real: low-order part of UT1 Julian date\ndelta_t::Real: Difference TT-UT1 at jd_high + jd_low in seconds\n\nOptional arguments\n\ngst_type: Return results as mean (:mean) or apparent (:apparent) time\nmethod: Computation method, CIO-based (:CIO) or equinox-based (:equinox)\naccuracy: Either :full or :reduced accuracy\n\n\n\n\n\n","category":"method"},{"location":"#NOVAS.tdb2tt-Tuple{Real}","page":"Home","title":"NOVAS.tdb2tt","text":"tdb2tt(tdb_jd)\n\nComputes the Terrestrial Time (TT) or Terrestrial Dynamic Time (TDT) Julian date corresponding to a Barrycentric Dynamical Time (TDB) Julian Date)\n\nArguments\n\ntdb_jd::Real: TDB Julian date\n\n\n\n\n\n","category":"method"}]
}
