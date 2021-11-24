using BenchmarkTools, DataFrames
import NOVAS
# Include un-exported c wrappers
include("../../test/wrapper.jl")

results = DataFrame(fname = [], c_time = [], julia_time = [])

# Utilities
c_fund_args = mean(@benchmark fund_args($rand())).time
julia_fund_args = mean(@benchmark NOVAS.fund_args($rand())).time
push!(results, ("fund_args", c_fund_args, julia_fund_args))

c_norm_ang = mean(@benchmark norm_ang($rand())).time
julia_norm_ang = mean(@benchmark NOVAS.norm_ang($rand())).time
push!(results, ("norm_ang", c_norm_ang, julia_norm_ang))

c_ee_ct_full = mean(@benchmark ee_ct($rand(), $rand(), 0)).time
julia_ee_ct_full = mean(@benchmark NOVAS.ee_ct($rand(), $rand(); accuracy = :full)).time
push!(results, ("ee_ct (full accuracy)", c_ee_ct_full, julia_ee_ct_full))

c_ee_ct_reduced = mean(@benchmark ee_ct($rand(), $rand(), 1)).time
julia_ee_ct_reduced = mean(@benchmark NOVAS.ee_ct($rand(), $rand(); accuracy = :reduced)).time
push!(results, ("ee_ct (reduced accuracy)", c_ee_ct_reduced, julia_ee_ct_reduced))

# Nutation
c_nu2000k = mean(@benchmark nu2000k($rand(), $rand())).time
julia_nu2000k = mean(@benchmark NOVAS.nu2000k($rand(), $rand())).time
push!(results, ("nu2000k", c_nu2000k, julia_nu2000k))

c_iau2000a = mean(@benchmark iau2000a($rand(), $rand())).time
julia_iau2000a = mean(@benchmark NOVAS.iau2000a($rand(), $rand())).time
push!(results, ("iau2000a", c_iau2000a, julia_iau2000a))

# NOVAS
c_nutation_angles = mean(@benchmark nutation_angles($rand(), 0)).time
julia_nutation_angles = mean(@benchmark NOVAS.nutation_angles($rand(), accuracy = :full)).time
push!(results, ("nutation_angles", c_nutation_angles, julia_nutation_angles))

c_mean_obliq = mean(@benchmark mean_obliq($rand())).time
julia_mean_obliq = mean(@benchmark NOVAS.mean_obliq($rand())).time
push!(results, ("mean_obliq", c_mean_obliq, julia_mean_obliq))

c_frame_tie = mean(@benchmark frame_tie($rand(3), 0)).time
julia_frame_tie = mean(@benchmark NOVAS.frame_tie($rand(3), :icrs2dynamic)).time
push!(results, ("frame_tie", c_frame_tie, julia_frame_tie))

c_nutation = mean(@benchmark nutation($rand(), 0, 0, $rand(3))).time
julia_nutation = mean(@benchmark NOVAS.nutation($rand(), $rand(3); accuracy = :full, direction = :mean2true)).time
push!(results, ("nutation", c_nutation, julia_nutation))

c_precession = mean(@benchmark precession($rand(), $rand(3), NOVAS.T0)).time
julia_precession = mean(@benchmark NOVAS.precession($rand(), $rand(3), NOVAS.T0)).time
push!(results, ("precession", c_precession, julia_precession))

c_e_tilt = mean(@benchmark e_tilt($rand(), 0)).time
julia_e_tilt = mean(@benchmark NOVAS.e_tilt($rand(); accuracy = :full)).time
push!(results, ("e_tilt", c_e_tilt, julia_e_tilt))

c_ira_equinox_mean = mean(@benchmark ira_equinox($rand(), 0, 0)).time
julia_ira_equinox_mean = mean(@benchmark NOVAS.ira_equinox($rand(); accuracy = :full, equinox = :mean)).time
push!(results, ("ira_equinox (mean)", c_ira_equinox_mean, julia_ira_equinox_mean))

c_ira_equinox_true = mean(@benchmark ira_equinox($rand(), 1, 0)).time
julia_ira_equinox_true = mean(@benchmark NOVAS.ira_equinox($rand(); accuracy = :full, equinox = :true)).time
push!(results, ("ira_equinox (true)", c_ira_equinox_true, julia_ira_equinox_true))

c_cio_location = mean(@benchmark cio_location($rand(), 0)).time
julia_cio_location = mean(@benchmark NOVAS.cio_location($rand(); accuracy = :full)).time
push!(results, ("cio_location", c_cio_location, julia_cio_location))

c_cio_basis = mean(@benchmark cio_basis($rand(), $rand(), 2, 0)).time
julia_cio_basis = mean(@benchmark NOVAS.cio_basis($rand(), $rand(); accuracy = :full)).time
push!(results, ("cio_basis", c_cio_basis, julia_cio_basis))

c_era = mean(@benchmark era($rand(), $rand())).time
julia_era = mean(@benchmark NOVAS.era($rand(), $rand())).time
push!(results, ("era", c_era, julia_era))

c_tdb2tt = mean(@benchmark tdb2tt($rand())).time
julia_tdb2tt = mean(@benchmark NOVAS.tdb2tt($rand())).time
push!(results, ("tdb2tt", c_tdb2tt, julia_tdb2tt))

c_sidereal_time_mean_cio = mean(@benchmark sidereal_time($rand(), $rand(), $rand(), 0, 0, 0)).time
julia_sidereal_time_mean_cio = mean(@benchmark NOVAS.sidereal_time($rand(), $rand(), $rand(); gst_type = :mean, method = :CIO, accuracy = :full)).time
push!(results, ("sidereal_time (CIO) (mean)", c_sidereal_time_mean_cio, julia_sidereal_time_mean_cio))

c_sidereal_time_apparent_cio = mean(@benchmark sidereal_time($rand(), $rand(), $rand(), 1, 0, 0)).time
julia_sidereal_time_apparent_cio = mean(@benchmark NOVAS.sidereal_time($rand(), $rand(), $rand(); gst_type = :apparent, method = :CIO, accuracy = :full)).time
push!(results, ("sidereal_time (CIO) (apparent)", c_sidereal_time_apparent_cio, julia_sidereal_time_apparent_cio))

c_sidereal_time_mean_equinox = mean(@benchmark sidereal_time($rand(), $rand(), $rand(), 0, 1, 0)).time
julia_sidereal_time_mean_equinox = mean(@benchmark NOVAS.sidereal_time($rand(), $rand(), $rand(); gst_type = :mean, method = :equinox, accuracy = :full)).time
push!(results, ("sidereal_time (equinox) (mean)", c_sidereal_time_mean_equinox, julia_sidereal_time_mean_equinox))

c_sidereal_time_apparent_equinox = mean(@benchmark sidereal_time($rand(), $rand(), $rand(), 1, 1, 0)).time
julia_sidereal_time_apparent_equinox = mean(@benchmark NOVAS.sidereal_time($rand(), $rand(), $rand(); gst_type = :apparent, method = :equinox, accuracy = :full)).time
push!(results, ("sidereal_time (equinox) (apparent)", c_sidereal_time_apparent_equinox, julia_sidereal_time_apparent_equinox))