# NOVAS.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kiranshila.github.io/NOVAS.jl/dev)
[![Build Status](https://github.com/kiranshila/NOVAS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kiranshila/NOVAS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kiranshila/NOVAS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kiranshila/NOVAS.jl)

## What is this
A reimplementation of the Naval Observatory Vector Astrometry Software (NOVAS) in pure Julia.

## Why does this exist
So far, all of the heavy-lifting in astrometrics calculations in Julia end up calling back to C libraries ([erfa](https://github.com/JuliaAstro/ERFA.jl)) or are feature-incomplete. As we want to build up an extensive library of high-performance astronomy and astrophysics software in Julia, having the foundation in pure Julia is quite important. This would then allow JuliaAstro software to be completley composable with other libraries (like [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl)) and be differentiable with any of the world-class AD tools available.

## Why not rewrite SOFA/ERFA/SLALIB/LIBNOVA....
ERFA, the Essential Routines for Fundamental Astronomy, replicated SOFA's functionality with bug fixes and an updated BSD 3-clause license. The other libraries are either feature-incomplete or unmaintained.

Between SOFA and NOVAS, the libraries tend to match each other in percission and speed, with NOVAS [performing much faster in some cases](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/10707/107071Z/A-comparison-of-SOFA-and-NOVAS-astrometric-software-libraries/10.1117/12.2311800.full?SSO=1). As such, NOVAS seemed to be the better candidate to rewrite.

## Implementation Details
In most cases, the algorithms are copied verbatim from the C implementation. However, linear algebra is used when appropriate along with Julia-optimized functions (`rem2pi` versus `fmod(x,2*PI)`). Additionally, the entire library is written with a functional flavor in mind. The vast majority of NOVAS accepts pointers to update instead of returning values. For these functions, tuples are returned instead.

## Performance and Testing
All functions are tested against their counterparts in C using an unexported C wrapper library that utilizes [NOVAS_jll](https://github.com/JuliaBinaryWrappers/NOVAS_jll.jl). This binary artifact was compiled with `-O2` flags. This Julia implementation usually exceeds the performance of the C version by at least two-fold and sometimes much more.

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