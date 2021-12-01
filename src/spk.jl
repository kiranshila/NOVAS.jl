export state, register_spk!

import Base.convert

abstract type AbstractSPKType end

# Type 1: Modified Difference arrays
#FIXME
struct SPKMDA <: AbstractSPKType
    coeffs::Matrix{Float64}
    epochs::Vector{Float64}
end

function convert(::Type{SPKMDA}, data::Vector{Float64})
    N = Int(data[end])
    diff_line_coeffs = hcat([data[i:(i + 70)] for i in 1:71:(71N)]...)
    epochs = [data[i] for i in (71(N + 1)):((71(N + 1)) + N)]
    return SPKMDA(diff_line_coeffs, epochs)
end

function (r::SPKMDA)(x::Real)
    return pos, vel
end

function state(entry::SPKMDA, time::Real)
    # Find the coeff line that is in range and evaluate
end

# Type 2: Chebyshev (position only)
struct ChebyshevRecord
    midpoint::Float64
    radius::Float64
    x::ChebyshevT
    y::ChebyshevT
    z::ChebyshevT
    ẋ::ChebyshevT
    ẏ::ChebyshevT
    ż::ChebyshevT
end

function (r::ChebyshevRecord)(x::Real)
    s = (x - r.midpoint) / r.radius
    pos = [r.x(s), r.y(s), r.z(s)]
    vel = [r.ẋ(s), r.ẏ(s), r.ż(s)]
    return pos, vel
end
struct SPKChebyshev <: AbstractSPKType
    records::Vector{ChebyshevRecord}
    INIT::Float64
    INTLEN::Float64
    RSIZE::Float64
    N::Float64
end

function convert(::Type{SPKChebyshev}, data::Vector{Float64})
    init = Int(data[end - 3])
    intlen = Int(data[end - 2])
    rsize = Int(data[end - 1])
    n = Int(data[end])
    degree = (rsize - 2) ÷ 3 - 1
    records = ChebyshevRecord[]
    @inbounds for record in 1:n
        offset = (record - 1) * rsize + 1
        mid = data[offset]
        radius = data[offset + 1]
        # Grab coefficients
        x_coefs = data[(offset + 2):(offset + 2 + degree)]
        y_coefs = data[(offset + 3 + degree):(offset + 3 + 2degree)]
        z_coefs = data[(offset + 4 + 2degree):(offset + 4 + 3degree)]
        # Build polynomials
        x = ChebyshevT(x_coefs)
        y = ChebyshevT(y_coefs)
        z = ChebyshevT(z_coefs)
        ẋ = derivative(x)
        ẏ = derivative(y)
        ż = derivative(x)
        # And create
        push!(records, ChebyshevRecord(mid, radius, x, y, z, ẋ, ẏ, ż))
    end
    # Sort
    sort!(records; by=x -> x.midpoint)
    return SPKChebyshev(records, init, intlen, rsize, n)
end

function state(entry::SPKChebyshev, time::Real)
    # Verify the target time is even in this range
    ending = entry.INIT + entry.INTLEN * entry.N
    @assert entry.INIT <= time <= ending "Target time not in range"
    # Binary search through segments to find the polynomial we care about
    found = false
    search_range = 1:Int(entry.N)
    current = -1
    while !found
        current = (search_range[end] - search_range[1]) ÷ 2 + search_range[1]
        # Check if time is covered by this segment
        current_record = entry.records[current]
        rstart = current_record.midpoint - current_record.radius
        rend = current_record.midpoint + current_record.radius
        if rstart <= time <= rend
            found = true
            break
        elseif length(search_range) == 1
            break
        elseif time < rstart
            search_range = search_range[1]:(current - 1)
        elseif time > rend
            search_range = (current + 1):search_range[end]
        end
    end
    @assert found "Segments don't cover requested time"
    # Evaluate time on segment
    return entry.records[current](time)
end

# There are more types, do we need them?
const SPK_TYPES = Dict(1 => SPKMDA, 2 => SPKChebyshev)
const SPK_REGISTRY = Dict{NTuple{3,Int},AbstractSPKType}()

"""
    register_spk!(daf)

Registers the ephemeris data from a DAF file into the global registry
"""
function register_spk!(daf::DAFFile)
    # Reinterpret DAF into an SPK file
    @assert daf.id == "DAF/SPK " "Input DAF is not a valid SPK file"
    # Each array in the DAF file represents a segment for a given descriptor
    # We are going to represent it as (target,center,reference) => callable
    for segment in daf.arrays
        target = segment.ints[1]
        center = segment.ints[2]
        reference = segment.ints[3]
        @assert segment.ints[4] ∈ keys(SPK_TYPES) "Unknown SPK representation"
        type = SPK_TYPES[segment.ints[4]]
        callable = convert(type, segment.data)
        SPK_REGISTRY[(target, center, reference)] = callable
    end
end

function register_spk!(filename::String)
    return register_spk!(read_daf(filename))
end

"""
    state(time,target,center)

Evaluates the state vectors (position and velocity) at a given time for a target and center system.

# Arguments
- `time::Real`: Time in TDB seconds since J2000
- `target::Int`: NAIF code for target
- `center::Int`: NAIF code for center

# Optional Arguments
- `reference::Int`: NAIF code for reference system
"""
function state(time::Real, target::Int, center::Int, reference::Int=1)
    # Check to see if a record is registered with this combination
    key = (target, center, reference)
    @assert key ∈ keys(SPK_REGISTRY) "No SPK entry found for this target/center/reference combination"
    # Evaluate state for the one that matches
    return state(SPK_REGISTRY[key], time)
end

function state(time::Real, target::String, center::String, reference::String="J2000")
    return state(time, naif_body_codes[target], naif_body_codes[center],
                 naif_reference_codes[reference])
end
