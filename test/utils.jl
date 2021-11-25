using RandomizedPropertyTest, NOVAS, BenchmarkTools, Random, DataFrames

rng = MersenneTwister(0)

# Include un-exported c wrappers
include("wrapper.jl")

# Add custom isapprox to test tuples
import Base.isapprox
Base.isapprox(x::Tuple, y::Tuple; kws...) = isapprox(collect(x), collect(y); kws...)

# Setup record of benchmark results
results = DataFrame(fname = String[], c_time = Float64[], julia_time = Float64[], speedup = Float64[])

# Benchmark using curried f from RPT
function quickbench(f::Function, types)
    genargs = () -> RandomizedPropertyTest.generate(rng, types)
    median(@benchmark $f($genargs()...)).time
end

# Stolen from RandomizedPropertyTest
macro testbench(args...)
    names = Symbol[]
    types = []
    nexpr = 10^4
    expr = args[1]
    vartypes = args[2:length(args)]
    for e in vartypes
        if e.args[1] isa Symbol
            newsymbols = Symbol[e.args[1]]
        elseif e.args[1].head == :tuple && all(x -> x isa Symbol, e.args[1].args)
            newsymbols = e.args[1].args
        end
        for symb in newsymbols
            push!(names, symb)
            push!(types, e.args[2])
        end
    end
    nametuple = Expr(:tuple, names...)
    typetuple = esc(Expr(:tuple, types...))
    exprstr = let io = IOBuffer()
        print(io, expr)
        seek(io, 0)
        read(io, String)
    end
    namestrs = [String(n) for n in names]
    # Build equivalence test
    NOVAS_expr = :(NOVAS.$(expr.args[1])($(expr.args[2:end]...)))
    test_expr = :($expr ≈ $NOVAS_expr)
    test_fexpr = esc(Expr(:(->), nametuple, test_expr))
    c_fexpr = esc(Expr(:(->), nametuple, expr))
    j_fexpr = esc(Expr(:(->), nametuple, NOVAS_expr))
    return quote
        test_result = RandomizedPropertyTest.do_quickcheck($test_fexpr, $exprstr, $namestrs, $typetuple, $nexpr)
        if test_result
            ctime = quickbench($c_fexpr, $typetuple)
            jtime = quickbench($j_fexpr, $typetuple)
            speedup = ctime / jtime
            push!(results, ($exprstr, ctime, jtime, speedup))
        end
        test_result
    end
end

# Custom types for common arguments
struct Position end
RandomizedPropertyTest.generate(rng::AbstractRNG, _::Type{Position}) = rand(Float64, 3)

struct Location end

function RandomizedPropertyTest.generate(rng::AbstractRNG, _::Type{Location})
    lat = rand(Cdouble) * 180 - 90
    lon = rand(Cdouble) * 180 - 90
    alt = rand(Cdouble) * 8000
    temp = rand(Cdouble) * 50
    pressure = rand(Cdouble) * 1000
    NOVAS.OnSurface(lat, lon, alt, temp, pressure)
end