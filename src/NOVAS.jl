module NOVAS

using DataDeps, LinearAlgebra, Memoize, StaticArrays, DelimitedFiles
import Base.show

include("init.jl")
include("constants.jl")
include("utils.jl")
include("nutation.jl")
include("solarsystem.jl")
include("main.jl")
# Thin C wrapper
include("LibNOVAS.jl")

end
