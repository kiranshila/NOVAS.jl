# Register and Download DataDeps *before* the rest of modulue init, as we need them for the const declarations

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
    post_fetch_method = unpack
))
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
    """,
    "https://github.com/kiranshila/nu2000k/archive/main.zip",
    "48b070d61c464da2abcea6a2221831e8ddd78b8d611bb4d9bcc712a50be9eeff";
    post_fetch_method = function(fn)
        unpack(fn)
        dir = "nu2000k-main"
        innerfiles=  readdir(dir)
        mv.(joinpath.(dir,innerfiles),innerfiles)
        rm(dir,recursive=true)
        rm("main.zip")
    end))
# Eagerly grab data deps
datadep"2003IERSConventions"
datadep"nu2000k"