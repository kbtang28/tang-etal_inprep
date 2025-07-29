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

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))

pt = 4/3
inch = 96

# plot A-B, B-C, and I-J IF utilization for limited wind & solar resource scenario
fig = Figure(size=(0inch, 3.25inch), fontsize=10pt)
ls_df = scale_load_shed_results(69, 2011)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls)

if_flow_df = scale_if_flow_results(69, 2011)

ax1 = Axis(fig[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)),
           xlabel="hour", ylabel="utilization",
           title="(a) A-B interface utilization",
           titlegap=5.0
)
hidexdecorations!(ax1)
vspan!(ax1, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax1, if_flow_df[1, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

ax2 = Axis(fig[2,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)),
           xlabel="hour", ylabel="utilization",
           title="(b) B-C interface utilization",
           titlegap=5.0
)
hidexdecorations!(ax2)
vspan!(ax2, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax2, if_flow_df[2, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

ax3 = Axis(fig[3,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)),
           xlabel="hour", ylabel="utilization",
           title="(c) I-J interface utilization",
           titlegap=5.0
)
vspan!(ax3, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax3, if_flow_df[10, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

linkxaxes!(ax1, ax3)
linkxaxes!(ax2, ax3)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace
ax3.yticklabelspace = yspace

colsize!(fig.layout, 1, Aspect(1, 12.0))
rowgap!(fig.layout, 10.0)
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "s69_ABif_BCif_IJif_utilization.png"), fig)

# plot zone J wind availability vs. I-J IF utilization for limited wind & solar resource scenario
fig = Figure(size=(0inch, 4inch), fontsize=10pt)
_, _, _, wind_cap_sf, solar_cap_sf, _, du_scenario = get_du_factors(69)
ls_df = scale_load_shed_results(69, 2011)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
ls01 = .!isapprox.(ls, 0.0; atol=1e-2)

if_flow_df = scale_if_flow_results(69, 2011)
if10 = if_flow_df[10, :scaled_if_flow][JJA_idx]

flow_df = scale_flow_results(69,2011)
flow75 = flow_df[75, :scaled_flow][JJA_idx]

wind_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "wind", "wind_2011.csv"), ',', Float64)[:, 2:end]
wind_cap_mat = wind_cap_sf*round.(wind_cap_mat[:, JJA_idx]; digits=2)
wind_cap_J = wind_cap_mat[16, :]
wind_cap_J_theoretical = wind_cap_sf*8250
wind_cap_J_ratio = wind_cap_J./wind_cap_J_theoretical
axmain1 = Axis(fig[2,1], 
               xgridvisible=false, ygridvisible=false, 
               xlabel="I-J interface utilization", 
               ylabel="Prop. of zone J wind available",
               yticks=0.0:0.5:1.0, 
               limits=((0.0, 1.0), (-0.2, 1.0))
)
axmain2 = Axis(fig[2,2],
               xgridvisible=false, ygridvisible=false,
               xlabel="I-J branch utilization",
               ylabel="Prop. of zone J wind available",
               yticks=0.0:0.5:1.0,
               limits=((0.0, 1.0), (-0.2, 1.0))
)
gtop = fig[1,1:2] = GridLayout()
axtop1 = Axis(gtop[1,1])
axtop2 = Axis(gtop[1,2])
axright = Axis(fig[2,3])
linkxaxes!(axmain1, axtop1)
linkxaxes!(axmain2, axtop2)
linkyaxes!(axmain1, axright)
linkyaxes!(axmain2, axright)
hideydecorations!(axmain2)

vlines!(axmain1, [0.0], linestyle=:dash, color=:black)
vlines!(axmain2, [0.0], linestyle=:dash, color=:black)
vlines!(axtop1, [0.0], linestyle=:dash, color=:black)
vlines!(axtop2, [0.0], linestyle=:dash, color=:black)
for (idx, color, alpha) in zip([.!ls01, ls01], ["#7f7f7f", "#d32f2f"], [0.5, 0.3])
    scatter!(axmain1, 
             if10[idx], 
             wind_cap_J_ratio[idx], 
             color=(color, alpha), 
    )
    scatter!(axmain2, 
             flow75[idx], 
             wind_cap_J_ratio[idx], 
             color=(color, alpha), 
    )
    violin!(axtop1, 
            zeros(sum(idx)), if10[idx], 
            color=(color, alpha), 
            side=:right, orientation=:horizontal, 
            datalimits=extrema
    )
    violin!(axtop2,
            zeros(sum(idx)), flow75[idx],
            color=(color, alpha),
            side=:right, orientation=:horizontal,
            datalimits=extrema
    )
    violin!(axright, 
            zeros(sum(idx)), wind_cap_J_ratio[idx], 
            color=(color, alpha), 
            side=:right, 
            datalimits=extrema
    )
end
text!(axmain1, [maximum(if10)], [maximum(wind_cap_J_ratio)], text="(a)", font=:bold, align=(:center, :center))
text!(axmain2, [maximum(flow75)], [maximum(wind_cap_J_ratio)], text="(b)", font=:bold, align=(:center, :center))
hidedecorations!(axtop1, grid=true)
hidedecorations!(axtop2, grid=true)
hidedecorations!(axright, grid=true)

elem1 = PolyElement(color=("#d32f2f", 0.3), points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
elem2 = PolyElement(color=("#7f7f7f", 0.5), points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
leg = Legend(fig[1,3], 
             [elem1, elem2], ["zone G\nshortage", "no zone G\nshortage"],
             framevisible=false,
             padding=(13.0, 0.0, 0.0, 0.0),
             rowgap=10,
             patchlabelgap=8
)

rowsize!(fig.layout, 1, Relative(1/4))
rowsize!(gtop, 1, Aspect(1, 0.35))
colsize!(fig.layout, 1, Aspect(2, 1.0))
colsize!(fig.layout, 2, Aspect(2, 1.0))
colsize!(fig.layout, 3, Aspect(2, 0.35))
rowgap!(fig.layout, 1, 7)
colgap!(fig.layout, 5)
colgap!(gtop, 1, 5)
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "s69_zoneJwind_IJif_scatter.png"), fig, dpi=300)