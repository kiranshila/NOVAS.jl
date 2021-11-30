# This is the interface to the JPL ephemerides for use with eph_manager
# The origin of the source here is `solsys1.c`, as that is what is bundled with pyNOVAS


function planet_ephemeris(tjd::AbstractVector,target,center)

end

"""
    solarsystem(tjd,body)

Provides an interface between the JPL direct-access solar system ephemerides and NOVAS

# Arguments
- `tjd::Real`:Julian date of the desired time, on the TDB time scale
- `body::Int`: Body identification number

# Optional Arguments
- `origin::Symbol=:barycenter`: The originof calculations, solar system `:barycenter`, the center of mass of the `:sun`, or the center of the `:earth`.

# Returns
`(pos,vel)` where
- `pos`: Position vector of `body` at `tjd` in AU
- `vel`: Velocity vector of `body` at `tdj` in AU/day
"""
function solarsystem(tjd::Real,body::Int;
    origin::Symbol=:barycenter)
    # Perform sanity checks
    @assert 1 <= body <= 11 "Illegal body id"
    @assert origin ∈ Set([:barycenter, :sun, :earth])

end
