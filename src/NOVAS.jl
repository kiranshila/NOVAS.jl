module NOVAS

using DataDeps,
    LinearAlgebra,
    Memoize,
    StaticArrays,
    DelimitedFiles

include("init.jl")
include("constants.jl")
include("utils.jl")
include("nutation.jl")
include("novas.jl")
# Thin C wrapper
include("LibNOVAS.jl")

end