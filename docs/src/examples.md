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

Next, we will use the current Julian UT1 date.

```@example sidereal
ut1_now = from_utc(now(); scale=UT1)
```

Then, we will use `AstroTime` to convert this date to a Julian date and give us the current offset to Terrestrial Time

```@example sidereal
jd_high, jd_low = value.(julian_twopart(ut1_now))
ΔT = getoffset(ut1_now, TT)
```

Finally, we will call `sidereal_time` to find the Greenwich mean sidereal time using the CIO method with full accuracy

```@example sidereal
sidereal_time(jd_high, jd_low, ΔT)
```

And because this is Julia, we can easily compose with other libraries, like [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl).
For example, we can add [uncertainty](https://eclipse.gsfc.nasa.gov/SEcat5/uncertainty.html) data to ΔT and see how it impacts our results

```@example sidereal
using Measurements
sidereal_time(jd_high, jd_low, ΔT ± 3)
```

## Position of a Star

This is adapted from the NOVASC manual, Example 1 in Chapter 3.

Suppose we have the catalog data from star HD 103095 (Groombridge 1830) for epoch J2000.0, expressed in the ICRS. We find this by pulling up its entry in something like [SIMBAD](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Gmb%201830). We can pull in a library like [AstroAngles](https://github.com/JuliaAstro/AstroAngles.jl) to simplify some of the angle conversions.

```@setup star
using NOVAS
```

```@example star
using AstroAngles

ra = hms"11:52:58.7683801554"ha
dec = dms"37:43:7.240082865"deg
μα = 4002.567
μδ = -5817.856
plax = 108.9551
v_rad = -98.008
```

Unlike the C version of NOVAS, we don't need to call `make_cat_entry` as we can make an instance of the struct at runtime.

```@example star
gmb_1830 = CatEntry("GMB 1830", "HD", 103095, ra, dec, μα, μδ, plax, v_rad)
```

Then, we call `app_star` to compute the apparent place of the star at a date `jd_tt`.

## Solar System Diagram

```@setup orbits
using NOVAS
```

In this example, we will play with the ephemeris data from JPL by plotting some trajectories.

First, we will set up the plotting system

```@example orbits
using Plots, PlotThemes
plotlyjs()
plot3d(; grid=false, background_color=:black, axis=false);
```

Next, we have to register which ephemeris dataset we will use.

```@example orbits
register_spk!(de440())
```

Then, we will use NOVAS to get the positions of each planet, with respects to the solar system barycenter

```@example orbits
days = 60 * 60 * 24
year = days * 365
years = 0:days:(2year)

get_planet_position(planet) = hcat((first.(state.(years, planet, 0)))...)

planets = ["Mercury", "Venus", "Earth", "Mars"]

# Plot the Sun!
scatter!([0], [0], [0]; label="Solar System Barycenter", markersize=3, markershape=:diamond,
         markercolor=:yellow)

for (i, planet) in enumerate(planets)
    position = get_planet_position(i)
    plot3d!(position[1, :], position[2, :], position[3, :]; label=planet)
end
```

It might be cool to see our moon's trajectory as well, but its record is stored as a vector from the Earth-Moon barycenter to the Moon. We can just add these tip to tail to get the Moon's position relative to the solar system's barycenter.

```@example orbits
# Plot the moon
moon_from_earth = hcat((first.(state.(years, 301, 3)))...)
earth = get_planet_position(3)
moon = moon_from_earth + earth
plot3d!(moon[1, :], moon[2, :], moon[3, :]; label="Luna")
```
