# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# load dependencies and helper functions
using CairoMakie
using CSV
using DataFrames

import ColorSchemes: seaborn_colorblind

# read in data
scenario_features_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "output", "scenario_features.csv")))

# plot
pt = 4/3
inch = 96
fig = Figure(size=(5inch, 0inch), fontsize=10pt)

ax1 = Axis(fig[1,1], xgridvisible=false, ygridvisible=false, xlabel="temp. anomaly (°C)", ylabel="hydro avail. anomaly (10⁵ MWh)", xtickformat="{:.1f}", ytickformat="{:.1f}")
scatter!(ax1, scenario_features_df[:, :avg_temp_anomaly], scenario_features_df[:, :avg_hydro_anomaly]/1e5, alpha=0.3)
Label(fig[1,1, TopRight()], "(a)", font=:bold, halign=:left, padding=(-25, 0, -25, 0))

ax2 = Axis(fig[1,2], xgridvisible=false, ygridvisible=false, xlabel="temp. anomaly (°C)", ylabel="load anomaly (10⁴ MWh)", xtickformat="{:.1f}", ytickformat="{:.1f}")
scatter!(ax2, scenario_features_df[:, :avg_temp_anomaly], scenario_features_df[:, :avg_load_anomaly]/1e4, alpha=0.3)
Label(fig[1,2, TopLeft()], "(b)", font=:bold, halign=:right, padding=(0, -25, -25, 0))

rowsize!(fig.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "figures", "load_and_hydro_vs_temp.png"), fig, dpi=300)