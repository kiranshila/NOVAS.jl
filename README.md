# NOVAS.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kiranshila.github.io/NOVAS.jl/dev)
[![Build Status](https://github.com/kiranshila/NOVAS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kiranshila/NOVAS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kiranshila/NOVAS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kiranshila/NOVAS.jl)

## What is this

A reimplementation of the Naval Observatory Vector Astrometry Software (NOVAS) in pure Julia.

## Why does this exist

So far, all of the heavy-lifting in astrometrics calculations in Julia end up calling back to C libraries ([erfa](https://github.com/JuliaAstro/ERFA.jl)) or are feature-incomplete. As we want to build up an extensive library of high-performance astronomy and astrophysics software in Julia, having the foundation in pure Julia is quite important. This would then allow JuliaAstro software to be completley composable with other libraries (like [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl)) and be differentiable with any of the world-class AD tools available.

## License Info

The license situation with NOVAS is strange. As noted in the unusually hard to get C source, the NOVAS source does not come with *any* copyright.

From their source
> This software was produced by the United States Naval Observatory at the 
> expense of United States taxpayers, and is therefore not suseptible to
> copyright, because a copyright would place taxpayer property under
> private ownership. Since it is not copyrighted, it cannot be licensed;
> it is simply free.

That being the case, this software will be covered under the MIT License as the contributors are not developing this at the expense of the US Government, and to make any further licensing issues simpler. This seems within the terms of the original source license, considering this is indeed derivative work, but the usual IANAL caveats apply.

## Original NOVAS Authors

We would like to thank the NOVAS authors for their hard work in creating the excellent software for which this is based.

**Barron, E. G., Kaplan, G. H., Bangert, J., Bartlett, J. L., Puatua, W., Harris, W., & Barrett, P. (2011)** `"Naval Observatory Vector Astrometry Software (NOVAS) Version 3.1, Introducing a Python Edition," <http://aa.usno.navy.mil/software/novas/novas_py/novas.pdf>`_ **Bull. AAS, 43, 2011.**