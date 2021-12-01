export de440

register(DataDep("2003IERSConventions",
                 """
                 2003 IERS Conventions
                 Author: Dennis D. McCarthy and Gérard Petit.
                 Website: https://iers-conventions.obspm.fr/index.php

                 The IERS Conventions describes the reference systems realized by the IERS, in addition to 
                 developing and maintaining the models and procedures used to support this endeavor. 

                 Citation: IERS Conventions (2003). Dennis D. McCarthy and Gérard Petit. 
                         (IERS Technical Note ; 32) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 
                         2004. 127 pp., paperback, ISBN 3-89888-884-3 (print version)
                 """,
                 "https://iers-conventions.obspm.fr/archive/2003/2003IERSConventionsArchive.tar.gz",
                 "1c7ac5c46f49dbfc6e140a34228067b480796ac11e3348ab41905e5af119f708";
                 post_fetch_method=unpack))
# I would like to use a ManualDataDep here, but that doesn't work for unregistered packages for some reason
# Once we register, we will move those constants back here
register(DataDep("nu2000k",
                 """
                 NU2000K Nutation Model
                 Author: George H. Kaplan
                 Website: https://www.usno.navy.mil/USNO/astronomical-applications/software-products/novas/novas-c

                 The 488-term nu2000k nutation model utilizes the luni-solar fundamental argument, planetary longitude, and general 
                 precession in longitude expressions from Simon et al. (1994) throughout. The NU2000K Nutation Model is fully documented in the 
                 United States Naval Observatory's NOVAS-C 3.0 software package.
                 """, "https://github.com/kiranshila/nu2000k/archive/main.zip",
                 "48b070d61c464da2abcea6a2221831e8ddd78b8d611bb4d9bcc712a50be9eeff";
                 post_fetch_method=function (fn)
                     unpack(fn)
                     dir = "nu2000k-main"
                     innerfiles = readdir(dir)
                     mv.(joinpath.(dir, innerfiles), innerfiles)
                     rm(dir; recursive=true)
                     return rm("main.zip")
                 end))

"""
iers2003(file)

Retrives and gives the location of a `file` from the 2003 IERS Convention Chapter 5 datapack
"""
function iers2003(file::String)
    return resolve("2003IERSConventions/chapter5/$file", @__FILE__)
end

"""
nu2000k(file)

Retrives and gives the location of a `file` from the NU2000k nutation model
"""
function nu2000k(file::String)
    return resolve("nu2000k/$file", @__FILE__)
end

function __init__()
    return register(DataDep("de440",
                            """
                            The JPL Planetary and Lunar Ephemerides DE440 and DE441
                            Author: Ryan S. Park, William M. Folkner, James G. Williams, and Dale H. Boggs
                            Website: https://iopscience.iop.org/article/10.3847/1538-3881/abd414

                            The latest JPL ephemeris with fully consistent treatment of planetary and 
                            lunar laser ranging data is DE440 (Park et al., 2021). The dynamical model
                            for DE440 includes a frictional damping between the fluid core and the
                            elastic mantle. This damping term is not suitable for extrapolation more
                            than several centuries into the past. In order to cover a longer time span,
                            the ephemeris DE441 was integrated without the lunar core/mantle damping term.
                            The positions of the planets for DE441 agree with the positions on DE440 to 
                            within one meter over the time covered by DE440. For the Moon DE441 differs
                            from DE440 mainly in the estimated tidal damping term causing a difference
                            in along-track position of the Moon of ~10 meters 100 years from the present
                            and growing quadratically for times more than 100 years from present. 

                            Citation: Ryan S. Park et al 2021 AJ 161 105
                            """,
                            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
                            "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"))
end

"""
    de440()

Retrives and gives the location of the JPL DE440 Planetary Ephemeris
"""
function de440()
    return datadep"de440/de440.bsp"
end
