```@raw html
<script>
    require.config({
        paths: { plotly: 'https://cdn.plot.ly/plotly-2.6.3.min' },
        shim: { plotly: { exports: 'plotly' } }
    });
</script>
```

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

In this example, we will play around with visualizing some of the JPL planetary ephemeris data. To do so, and to have a nice interactive plot, we will pull in PlotlyJS.

```@setup orbits
import UUIDs
using NOVAS, PlotlyJS

function documenter_plotly(plt)
    uuid = UUIDs.uuid4()
    html = """
        <div id=\"$(uuid)\"></div>
        <script>
        require(['plotly'], function (plotly) {
            PLOT = document.getElementById('$(uuid)');
            console.log(plotly);
            plotly.newPlot(PLOT, $(string(plt.plot.data)), $(string(plt.plot.layout)), { scrollZoom: true });
        });
        </script>
    """
    HTML(html)
end
```

First, we need to register which ephemeris kernel we want to use with the global registry. We'll use the latest DE440 planetary kernel.

```@example orbits
register_spk!(de440())
```

Then, we will write a utility function to grab the position vector from the NOVAS [`state`](@ref) function as well as calculate the time span we want to look at, say two years from J2000. Notice we can just use broadcasting here over the range.

```@example orbits
days = 60 * 60 * 24
year = days * 365
years = 0:3days:(2year)

get_planet_position(planet) = hcat((first.(state.(years, planet, 0)))...)
```

Now, we can simply call that function for each planet we're interested in. Included is a small utility function to work with the PlotlyJS API.

```@example orbits
planets = ["Mercury", "Venus", "Earth", "Mars"]
traces = GenericTrace[]
trace_orbit(x, y, z, name) = scatter3d(; x=x, y=y, z=z, mode="lines", name=name)
for (i, planet) in enumerate(planets)
    position = get_planet_position(i)
    push!(traces, trace_orbit(position[1, :], position[2, :], position[3, :], planet))
end
```

To make things more interesting, we can plot the solar system barycenter, as a point of reference

```@example orbits
push!(traces,
      scatter3d(; x=[0], y=[0], z=[0], name="Solar System Barycenter",
                marker=attr(; color="yellow", size=5, symbol="diamond")))
nothing #hide
```

Now, we just setup our Plotly layout and look at the pretty curves!

```@example orbits
layout = Layout(; plot_bgcolor="black", paper_bgcolor="black",
                legend=attr(; font=attr(; color="white")), font=attr(; color="white"),
                scene=attr(; xaxis=attr(; visible=false), yaxis=attr(; visible=false),
                           zaxis=attr(; visible=false)))
plt = plot(traces, layout)
documenter_plotly(plt) #hide
```
