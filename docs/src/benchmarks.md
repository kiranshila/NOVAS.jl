# Benchmarks

Here is a comparison of the runtime speed of every function implemented.
This is run of GitHub CI with system specs of:

```@example
using InteractiveUtils
versioninfo()
```

## Results

Unless otherwise noted, benchmark times are with "full" accuracy.
This usually corresponds to using the iau2000a nutation model.

```@setup bench
include("../../test/benchmark.jl")
```


```@example bench
results
```