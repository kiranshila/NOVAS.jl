module NOVAS

using DataDeps, LinearAlgebra, Memoize, StaticArrays, DelimitedFiles, FortranFiles,
      Polynomials
import Base.show

include("data.jl")
include("constants.jl")
include("utils.jl")
include("nutation.jl")
include("daf.jl")
include("naif_codes.jl")
include("spk.jl")
include("solarsystem.jl")
include("main.jl")
# Thin C wrapper
include("LibNOVAS.jl")

end
