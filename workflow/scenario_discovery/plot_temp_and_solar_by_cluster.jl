# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

# load dependencies and helper functions
using CairoMakie
using CSV
using DataFrames
using DelimitedFiles

import ColorSchemes: seaborn_colorblind

# load data
scenario_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "output", "scenario_features.csv")))
labels = readdlm(joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "clustering", "cluster_labels_3.txt"), Int)
insertcols!(scenario_df, :label => vec(labels))

# plot
pt = 4/3
inch = 96
fig = Figure(size=(3.5inch, 0inch), fontsize=8pt)
cluster_colors = Dict(1 => (seaborn_colorblind[1], 1.0), 2 => (seaborn_colorblind[5], 1.0), 3 => (seaborn_colorblind[3], 1.0))
color = [cluster_colors[scenario_df.label[i]] for i in 1:nrow(scenario_df)]
ax1 = Axis(fig[1,2],
           xgridvisible=false, ygridvisible=false,
           xlabel="cluster", ylabel="solar availability (10³ MWh)",
           limits=((0.25, 3.75), nothing),
           xticks=1:3, ytickformat="{:.1f}"
)
boxplot!(ax1, 
         scenario_df.label, scenario_df.avg_solar_anomaly/1e3, 
         color=color
)

ax2 = Axis(fig[1,1],
           xgridvisible=false, ygridvisible=false,
           xlabel="cluster", ylabel="temperature change (°C)",
           limits=((0.25, 3.75), nothing),
           xticks=1:3, ytickformat="{:.1f}"
)
boxplot!(ax2, 
         scenario_df.label, scenario_df.avg_temp_anomaly, 
         color=color
)

rowsize!(fig.layout, 1, Aspect(1, 2.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "clusters_solar_temp_boxplots.png"), fig)