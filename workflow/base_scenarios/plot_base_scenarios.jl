# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

# load dependencies
using CairoMakie
using CSV
using DataFrames
using DelimitedFiles

import ColorSchemes: tab20, amp, seaborn_colorblind, seaborn_deep

include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))

# helper function
function get_ls_span(ls)
    # returns ranges for coloring LS hours with vspan
    n = length(ls)
    ls01 = .!isapprox.(ls, 0.0; atol=1e-2)
    xlow = Int[]
    xhigh = Int[]

    if ls01[1] == 1
        push!(xlow, 1)
    end

    for i in 2:n
        if ls01[i-1] == 0 && ls01[i] == 1
            push!(xlow, i)
        end

        if ls01[i-1] == 1 && ls01[i] == 0
            push!(xhigh, i)
        end
    end

    if length(xlow) == length(xhigh)
        return xlow, xhigh
    else
        push!(xhigh, n)
        return xlow, xhigh
    end
end

pt = 4/3
inch = 96

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))

# plot scenario features for full ensemble with base scenarios highlighted
fig = Figure(size=(0inch, 2.25inch), fontsize=8pt)
scenario_features_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "output", "scenario_features.csv")))
scen_yr_pairs = [(140, 2012), (69, 2011), (290, 2002)]
ids = [findall((scenario_features_df.scenario .== scen) .&& (scenario_features_df.yr .== yr))[1] for (scen, yr) in scen_yr_pairs]
ax1 = Axis(fig[1,1], 
           xgridvisible=false, ygridvisible=false, 
           xlabel="solar avail. anomaly (10³ MWh)", ylabel="wind avail. anomaly (10³ MWh)", 
           xtickformat="{:.1f}", ytickformat="{:.1f}",
           xticklabelsize=7pt, yticklabelsize=7pt
)
scatter!(ax1, 
         scenario_features_df.avg_solar_anomaly/1e3, scenario_features_df.avg_wind_anomaly/1e3,
         markersize=6.0, 
         color="#a1a1a1", alpha=0.2
)
scatter!(ax1, 
         scenario_features_df.avg_solar_anomaly[ids]/1e3, scenario_features_df.avg_wind_anomaly[ids]/1e3, 
         color=[seaborn_colorblind[3], "#5d55cd", seaborn_colorblind[2]], 
         strokecolor=:black, strokewidth=1.5, 
         markersize=10.0
)
text!(ax1, 
      [maximum(scenario_features_df.avg_solar_anomaly/1e3)], [maximum(scenario_features_df.avg_wind_anomaly/1e3)], 
      text="(a)", 
      font=:bold, align=(:right, :top)
)
text!(ax1,
      scenario_features_df.avg_solar_anomaly[ids[1]]/1e3, scenario_features_df.avg_wind_anomaly[ids[1]]/1e3,
      text=["1"],
      align=(:right, :center),
      offset=(-10.0, 0.0)
)
text!(ax1,
      scenario_features_df.avg_solar_anomaly[ids[2]]/1e3, scenario_features_df.avg_wind_anomaly[ids[2]]/1e3,
      text=["2"],
      align=(:left, :center),
      offset=(10.0, 0.0)
)
text!(ax1,
      scenario_features_df.avg_solar_anomaly[ids[3]]/1e3, scenario_features_df.avg_wind_anomaly[ids[3]]/1e3,
      text=["3"],
      align=(:right, :center),
      offset=(-10.0, 0.0)
)

ax2 = Axis(fig[1,2], 
           xgridvisible=false, ygridvisible=false, 
           xlabel="temp. anomaly (°C)", ylabel="hydro avail. anomaly (10⁵ MWh)", 
           xtickformat="{:.1f}", ytickformat="{:.1f}",
           xticklabelsize=7pt, yticklabelsize=7pt
)
scatter!(ax2, 
         scenario_features_df.avg_temp_anomaly, scenario_features_df.avg_hydro_anomaly/1e5,
         markersize=6.0, 
         color="#a1a1a1", alpha=0.2
)
scatter!(ax2, 
         scenario_features_df.avg_temp_anomaly[ids], scenario_features_df.avg_hydro_anomaly[ids]/1e5, 
         color=[seaborn_colorblind[3], "#5d55cd", seaborn_colorblind[2]], 
         strokecolor=:black, strokewidth=1.5, 
         markersize=10.0
)
text!(ax2, 
      [maximum(scenario_features_df.avg_temp_anomaly)], [maximum(scenario_features_df.avg_hydro_anomaly/1e5)], 
      text="(b)", 
      font=:bold, align=(:right, :top)
)
text!(ax2,
      scenario_features_df.avg_temp_anomaly[ids], scenario_features_df.avg_hydro_anomaly[ids]/1e5,
      text=["1", "2", "3"],
      align=(:right, :center),
      offset=(-10.0, 0.0)
)

ax3 = Axis(fig[1,3], 
           xgridvisible=false, ygridvisible=false, 
           xlabel="prop. hours with E-G IF congestion", ylabel="prop. wind and solar curtailed", 
           xtickformat="{:.1f}", ytickformat="{:.1f}",
           xticklabelsize=7pt, yticklabelsize=7pt
)
scatter!(ax3, 
         scenario_features_df.prop_congested_pos_hrs, scenario_features_df.prop_curtailed,
         markersize=6.0, 
         color="#a1a1a1", alpha=0.2
)
scatter!(ax3, 
         scenario_features_df.prop_congested_pos_hrs[ids], scenario_features_df.prop_curtailed[ids], 
         color=[seaborn_colorblind[3], "#5d55cd", seaborn_colorblind[2]], 
         strokecolor=:black, strokewidth=1.5, 
         markersize=10.0
)
text!(ax3, 
      [maximum(scenario_features_df.prop_congested_pos_hrs)], [maximum(scenario_features_df.prop_curtailed)], 
      text="(c)", 
      font=:bold, align=(:right, :top)
)
text!(ax3,
      scenario_features_df.prop_congested_pos_hrs[ids], scenario_features_df.prop_curtailed[ids],
      text=["1", "2", "3"],
      align=(:right, :center),
      offset=(-10.0, 0.0)
)

elems = [MarkerElement(color=seaborn_colorblind[3], marker=:circle, strokewidth=1.5, markersize=10.0), 
         MarkerElement(color="#5d55cd", marker=:circle, strokewidth=1.5, markersize=10.0), 
         MarkerElement(color=seaborn_colorblind[2], marker=:circle, strokewidth=1.5, markersize=10.0)
]
axislegend(ax3, 
           elems, ["Well-behaved", "Limited wind\nand solar", "Extreme temp."], 
           position=:lb, rowgap=-0.1, patchlabelgap=2.0, padding=(0.0, 6.0, 0.0, 0.0),
           fontsize=7pt
)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "scenario_features.png"), fig, px_per_unit=300/inch)

# plot renewable capacity ratios
fig = Figure(size=(0inch, 2.8inch), fontsize=8pt)
# baseline scenario (scenario 140, 2012)
_, _, _, base_wind_cap_sf, base_solar_cap_sf, _, base_du_scenario = get_du_factors(140)
base_ls_df = scale_load_shed_results(140, 2012)
base_ls = base_ls_df[44, :scaled_load_shed][JJA_idx]
base_ls01 = .!isapprox.(base_ls, 0.0; atol=1e-2)
base_side = map(base_ls01) do ls
    return ls == true ? :right : :left
end

base_wind_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "wind", "wind_2012.csv"), ',', Float64)[:, 2:end]
base_wind_cap_mat = base_wind_cap_sf*round.(base_wind_cap_mat[:, JJA_idx]; digits=2)
base_wind_cap_C = base_wind_cap_mat[7, :]
base_wind_cap_C_theoretical = 0.5*1923*base_wind_cap_sf
base_wind_ratio_C = base_wind_cap_C./base_wind_cap_C_theoretical
base_wind_cap_E = base_wind_cap_mat[10, :]
base_wind_cap_E_theoretical = 0.7*1821*base_wind_cap_sf
base_wind_ratio_E = base_wind_cap_E./base_wind_cap_E_theoretical
base_wind_cap_K = vec(sum(base_wind_cap_mat[17:18, :], dims=1))
base_wind_cap_K_theoretical = 121*base_wind_cap_sf
base_wind_ratio_K = base_wind_cap_K./base_wind_cap_K_theoretical

base_solar_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "solar", "Scenario$(base_du_scenario)", "solarUPV2012.csv"), ',', Float64)[:, 2:end]
base_solar_cap_mat = base_solar_cap_sf*round.(base_solar_cap_mat[:, JJA_idx]; digits=2)
base_solar_cap_A = vec(sum(base_solar_cap_mat[3:4, :], dims=1))
base_solar_cap_A_theoretical = 0.3*14440*base_solar_cap_sf
base_solar_ratio_A = base_solar_cap_A./base_solar_cap_A_theoretical
base_solar_cap_B = base_solar_cap_mat[5, :]
base_solar_cap_B_theoretical = 1648*base_solar_cap_sf
base_solar_ratio_B = base_solar_cap_B./base_solar_cap_B_theoretical

gab = fig[1,1] = GridLayout()
ga = gab[1,1] = GridLayout()
axa = Axis(ga[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=(nothing, (0.0, 1.0)),
           xticks=(1:5, ["Zone A\nsolar", "Zone B\nsolar", "Zone C\nwind", "Zone E\nwind", "Zone K\nwind"]),
           ylabel="prop. of installed\ncapacity available",
           yticklabelsize=7pt,
           title="(a) Well-behaved scenario",
           titlegap=6.0
)
for (i, ratio) in enumerate([base_solar_ratio_A, base_solar_ratio_B])
    color = map(base_side) do s
        return s == :left ? "#eec680" : seaborn_colorblind[2]
    end
    violin!(axa, i*ones(length(JJA_idx)), ratio, side=base_side, datalimits=(0.0, 1.0), scale=:count, color=color, strokewidth=0.75, strokecolor=:black)
end
for (i, ratio) in enumerate([base_wind_ratio_C, base_wind_ratio_E, base_wind_ratio_K])
    color = map(base_side) do s
        return s == :left ? "#99d8c7" : seaborn_colorblind[3]
    end
    violin!(axa, (i+2)*ones(length(JJA_idx)), ratio, side=base_side, datalimits=(0.0, 1.0), scale=:count, color=color, strokewidth=0.75, strokecolor=:black)
end
ylims!(axa, low=0, high=1.0)

# limited resource scenario (scenario 69, 2011)
_, _, _, lr_wind_cap_sf, lr_solar_cap_sf, _, lr_du_scenario = get_du_factors(69)
lr_ls_df = scale_load_shed_results(69, 2011)
lr_ls = lr_ls_df[44, :scaled_load_shed][JJA_idx]
lr_ls01 = .!isapprox.(lr_ls, 0.0; atol=1e-2)
lr_side = map(lr_ls01) do ls
    return ls == true ? :right : :left
end

lr_wind_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "wind", "wind_2011.csv"), ',', Float64)[:, 2:end]
lr_wind_cap_mat = lr_wind_cap_sf*round.(lr_wind_cap_mat[:, JJA_idx]; digits=2)
lr_wind_cap_A = lr_wind_cap_mat[1, :]
lr_wind_cap_A_theoretical = 0.24*2692*lr_wind_cap_sf
lr_wind_ratio_A = lr_wind_cap_A./lr_wind_cap_A_theoretical

lr_solar_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "solar", "Scenario$(lr_du_scenario)", "solarUPV2011.csv"), ',', Float64)[:, 2:end]
lr_solar_cap_mat = lr_solar_cap_sf*round.(lr_solar_cap_mat[:, JJA_idx]; digits=2)
lr_solar_cap_B = lr_solar_cap_mat[5, :]
lr_solar_cap_B_theoretical = 1648*lr_solar_cap_sf
lr_solar_ratio_B = lr_solar_cap_B./lr_solar_cap_B_theoretical

gb = gab[1,2] = GridLayout()
axb = Axis(gb[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=((0.5, 2.65), (0.0, 1.0)),
           xticks=(1:2, ["Zone B\nsolar", "Zone A\nwind"]),
           yticklabelsize=7pt,
           title="(b) Limited resource scenario",
           titlegap=6.0
)
for (i, ratio) in enumerate([lr_solar_ratio_B])
    color = map(lr_side) do s
        return s == :left ? "#eec680" : seaborn_colorblind[2]
    end
    violin!(axb, i*ones(length(JJA_idx)), ratio, side=lr_side, datalimits=(0.0, 1.0), scale=:count, color=color, strokewidth=0.75, strokecolor=:black)
end
for (i, ratio) in enumerate([lr_wind_ratio_A])
    color = map(lr_side) do s
        return s == :left ? "#99d8c7" : seaborn_colorblind[3]
    end
    violin!(axb, (i+1)*ones(length(JJA_idx)), ratio, side=lr_side, datalimits=(0.0, 1.0), scale=:count, color=color, strokewidth=0.75, strokecolor=:black)
end
ylims!(axb, low=0, high=1.0)
hideydecorations!(axb)

# extreme temperature scenario (scenario 290, 2002)
_, _, _, ext_wind_cap_sf, ext_solar_cap_sf, _, ext_du_scenario = get_du_factors(290)
ext_ls_df = scale_load_shed_results(290, 2002)
ext_ls = ext_ls_df[44, :scaled_load_shed][JJA_idx]
ext_ls01 = .!isapprox.(ext_ls, 0.0; atol=1e-2)
ext_side = map(ext_ls01) do ls
    return ls == true ? :right : :left
end

ext_wind_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "wind", "wind_2002.csv"), ',', Float64)[:, 2:end]
ext_wind_cap_mat = ext_wind_cap_sf*round.(ext_wind_cap_mat[:, JJA_idx]; digits=2)
ext_wind_cap_C = ext_wind_cap_mat[8, :]
ext_wind_cap_C_theoretical = 0.05*1923*ext_wind_cap_sf
ext_wind_ratio_C = ext_wind_cap_C./ext_wind_cap_C_theoretical
ext_wind_cap_D = ext_wind_cap_mat[9, :]
ext_wind_cap_D_theoretical = 1935*ext_wind_cap_sf
ext_wind_ratio_D = ext_wind_cap_D./ext_wind_cap_D_theoretical
ext_wind_cap_G = vec(sum(ext_wind_cap_mat[13:14, :], dims=1))
ext_wind_cap_G_theoretical = 606*ext_wind_cap_sf
ext_wind_ratio_G = ext_wind_cap_G./ext_wind_cap_G_theoretical
ext_wind_cap_H = ext_wind_cap_mat[15, :]
ext_wind_cap_H_theoretical = 303*ext_wind_cap_sf
ext_wind_ratio_H = ext_wind_cap_H./ext_wind_cap_H_theoretical
ext_wind_cap_K = vec(sum(ext_wind_cap_mat[17:19, :], dims=1))
ext_wind_cap_K_theoretical = (121 + 0.5*6488)*ext_wind_cap_sf
ext_wind_ratio_K = ext_wind_cap_K./ext_wind_cap_K_theoretical

ext_solar_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "solar", "Scenario$(ext_du_scenario)", "solarUPV2002.csv"), ',', Float64)[:, 2:end]
ext_solar_cap_mat = ext_solar_cap_sf*round.(ext_solar_cap_mat[:, JJA_idx]; digits=2)
ext_solar_cap_A = ext_solar_cap_mat[2, :]
ext_solar_cap_A_theoretical = 0.2*14440*ext_solar_cap_sf
ext_solar_ratio_A = ext_solar_cap_A./ext_solar_cap_A_theoretical
ext_solar_cap_B = ext_solar_cap_mat[5, :]
ext_solar_cap_B_theoretical = 1648*ext_solar_cap_sf
ext_solar_ratio_B = ext_solar_cap_B./ext_solar_cap_B_theoretical
ext_solar_cap_G = ext_solar_cap_mat[12, :]
ext_solar_cap_G_theoretical = 0.7*3353*ext_solar_cap_sf
ext_solar_ratio_G = ext_solar_cap_G./ext_solar_cap_G_theoretical
ext_solar_cap_K = ext_solar_cap_mat[14, :]
ext_solar_cap_K_theoretical = 0.5*1441*ext_solar_cap_sf
ext_solar_ratio_K = ext_solar_cap_K./ext_solar_cap_K_theoretical

gc = fig[2,1] = GridLayout()
axc = Axis(gc[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=(nothing, (0.0, 1.0)),
           xticks=(1:9, ["Zone A\nsolar", "Zone B\nsolar", "Zone G\nsolar", "Zone K\nsolar", "Zone C\nwind", "Zone D\nwind", "Zone G\nwind", "Zone H\nwind", "Zone K\nwind"]),
           ylabel="prop. of installed\ncapacity available",
           yticklabelsize=7pt,
           title="(c) Extreme temperature scenario",
           titlegap=6.0
)
for (i, ratio) in enumerate([ext_solar_ratio_A, ext_solar_ratio_B, ext_solar_ratio_G, ext_solar_ratio_K])
    color = map(ext_side) do s
        return s == :left ? "#eec680" : seaborn_colorblind[2]
    end
    violin!(axc, i*ones(length(JJA_idx)), ratio, side=ext_side, datalimits=(0.0, 1.0), scale=:count, color=color, strokewidth=0.75, strokecolor=:black)
end
for (i, ratio) in enumerate([ext_wind_ratio_C, ext_wind_ratio_D, ext_wind_ratio_G, ext_wind_ratio_H, ext_wind_ratio_K])
    color = map(ext_side) do s
        return s == :left ? "#99d8c7" : seaborn_colorblind[3]
    end
    violin!(axc, (i+4)*ones(length(JJA_idx)), ratio, side=ext_side, datalimits=(0.0, 1.0), scale=:count, color=color, strokewidth=0.75, strokecolor=:black)
end
ylims!(axc, low=0, high=1.0)

elem1 = [PolyElement(color="#eec680", strokecolor=:black, strokewidth=1.0, points = Point2f[(0, 0), (1, 0), (0, 1)]),
         PolyElement(color="#99d8c7", strokecolor=:black, strokewidth=1.0, points = Point2f[(1, 1), (1, 0), (0, 1)])
]
elem2 = [PolyElement(color=seaborn_colorblind[2], strokecolor=:black, strokewidth=1.0, points = Point2f[(0, 0), (1, 0), (0, 1)]),
         PolyElement(color=seaborn_colorblind[3], strokecolor=:black, strokewidth=1.0, points = Point2f[(1, 1), (1, 0), (0, 1)])
]
Legend(gab[1,3], 
       [elem1, elem2], ["without power\nshortage", "with power\nshortage"], 
       orientation=:vertical, 
       tellwidth=false, tellheight=false, 
       patchsize=(20.0, 20.0), 
       patchlabelgap=10.0, 
       rowgap=13,
       framevisible=false
)

colsize!(fig.layout, 1, Aspect(2, 8.0))
colsize!(gab, 1, Relative(5/9))
colsize!(gab, 2, Relative(2.5/9))
colsize!(gab, 3, Relative(1.5/9))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "renew_ratios.png"), fig)

# plot curtailment and power shortages for all three scenarios
fig = Figure(size=(10inch, 7inch), fontsize=10pt)
# well-behaved scenario
curtail_df = scale_curtail_results(140, 2012)
ls_df = scale_load_shed_results(140, 2012)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls)

ax1 = Axis(fig[1,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour",
           title="(a) Well-behaved scenario",
           titlegap=5.0
)
hidexdecorations!(ax1)
vspan!(ax1, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax1, curtail_df[7, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[1,3], color="#99d8c7")
Label(fig[1,3], "wind (C)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))
Label(fig[1:2,1], "curtailment", rotation=pi/2, tellheight=false)

ax2 = Axis(fig[2,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour"
)
vspan!(ax2, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax2, curtail_df[25, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[2,3], color="#eec680")
Label(fig[2,3], "solar (B)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))

# limited wind and solar resource scenario
curtail_df = scale_curtail_results(69, 2011)
ls_df = scale_load_shed_results(69, 2011)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls)

ax3 = Axis(fig[3,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour",
           title="(b) Limited wind & solar resource scenario",
           titlegap=5.0
)
hidexdecorations!(ax3)
vspan!(ax3, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax3, curtail_df[1, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[3,3], color="#99d8c7")
Label(fig[3,3], "wind (A)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))
Label(fig[3:4,1], "curtailment", rotation=pi/2, tellheight=false)

ax4 = Axis(fig[4,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour"
)
vspan!(ax4, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax4, curtail_df[25, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[4,3], color="#eec680")
Label(fig[4,3], "solar (B)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))

# extreme temperature scenario
curtail_df = scale_curtail_results(290, 2002)
ls_df = scale_load_shed_results(290, 2002)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls)

ax5 = Axis(fig[5,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour",
           title="(c) Extreme temperature scenario",
           titlegap=5.0
)
hidexdecorations!(ax5)
vspan!(ax5, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax5, curtail_df[18, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[5,3], color="#99d8c7")
Label(fig[5,3], "wind (J)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))
Label(fig[5:7,1], "curtailment", rotation=pi/2, tellheight=false)

ax6 = Axis(fig[6,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour"
)
hidexdecorations!(ax6)
vspan!(ax6, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax6, curtail_df[25, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[6,3], color="#eec680")
Label(fig[6,3], "solar (B)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))

ax7 = Axis(fig[7,2],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
           xlabel="hour"
)
vspan!(ax7, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax7, curtail_df[22, :scaled_curtail][JJA_idx], color=:gray22)
Box(fig[7,3], color="#eec680")
Label(fig[7,3], "solar (A)", rotation=pi/2, tellheight=false, fontsize=8pt, padding=(2.0, 2.0, 0, 0))

colgap!(fig.layout, 1, 5.0)
colgap!(fig.layout, 2, 0.0)
colsize!(fig.layout, 2, Aspect(1, 11.5))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "curtailment.png"), fig)