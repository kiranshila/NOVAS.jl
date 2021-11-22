# Benchmarks

Here is a comparison of the runtime speed of every function implemented.
This is run of GitHub CI with system specs of:

```@example
versioninfo()
```

## Results

```@setup bench
using BenchmarkTools, DataFrames
import NOVAS
# Include un-exported c wrappers
include("../test/wrapper.jl")

results = DataFrame(fname=[],c_time=[],julia_time=[])

# Utilities
c_fund_args      = mean(@benchmark fund_args($rand())).time
julia_fund_args  = mean(@benchmark NOVAS.fund_args($rand())).time
push!(results,("fund_args",c_fund_args,julia_fund_args))

c_norm_ang       = mean(@benchmark norm_ang($rand())).time
julia_norm_ang   = mean(@benchmark NOVAS.norm_ang($rand())).time
push!(results,("norm_ang",c_norm_ang,julia_norm_ang))

c_ee_ct_full     = mean(@benchmark ee_ct($rand(),$rand(),0)).time
julia_ee_ct_full = mean(@benchmark NOVAS.ee_ct($rand(),$rand();accuracy=:full)).time
push!(results,("ee_ct (full accuracy)",c_ee_ct_full,julia_ee_ct_full))

c_ee_ct_reduced     = mean(@benchmark ee_ct($rand(),$rand(),1)).time
julia_ee_ct_reduced = mean(@benchmark NOVAS.ee_ct($rand(),$rand();accuracy=:reduced)).time
push!(results,("ee_ct (reduced accuracy)",c_ee_ct_reduced,julia_ee_ct_reduced))

# Nutation
c_nu2000k     =  mean(@benchmark nu2000k($rand(),$rand())).time
julia_nu2000k =  mean(@benchmark NOVAS.nu2000k($rand(),$rand())).time 
push!(results,("nu2000k",c_nu2000k,julia_nu2000k))

c_iau2000a     =  mean(@benchmark iau2000a($rand(),$rand())).time
julia_iau2000a =  mean(@benchmark NOVAS.iau2000a($rand(),$rand())).time 
push!(results,("iau2000a",c_iau2000a,julia_iau2000a))
```

```@example bench
results
```