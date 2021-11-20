# Examples

## Calculating Sidereal Time
In this example, we will calculate the Greenwich mean sidereal time.

```@setup sidereal
using NOVAS
```
First, we need to pull in [AstroTime.jl](https://github.com/JuliaAstro/AstroTime.jl) to give us the most preceise time correction data and transformations.
We will also call it's update routine.
```@example sidereal
using AstroTime
AstroTime.update()
```
Next, we will use the current Julian UT1 date and find ΔT.
```@example sidereal
ut1_now = from_utc(now(); scale=UT1)
jd_high, jd_low = julian_twopart(ut1_now)
ΔT = getoffset(ut1_now,TT)
```
Finally, we will call `sidereal_time` to find the Greenwich mean sidereal time using the CIO method with full accuracy
```@example sidereal
sidereal_time(jd_high,jd_low,ΔT)
```