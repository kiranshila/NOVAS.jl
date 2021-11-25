# Benchmarks

Here is a comparison of the runtime speed of every function implemented.
This is run of GitHub CI with system specs of:

```@example
using InteractiveUtils
versioninfo()
```

## Results

```@setup bench
using CSV, DataFrames
```


```@example bench
println(pwd())
CSV.File("../../benchmarks.csv") |> DataFrame
```