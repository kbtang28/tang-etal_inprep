# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

# load dependencies and helper functions
using CairoMakie
using CSV
using DataFrames
using DelimitedFiles

import ColorSchemes: ColorScheme, seaborn_colorblind
import Colors: parse, Colorant
import StatsBase: mean

include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))

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
fig = Figure(size=(0inch, 3.25inch), fontsize=10pt)

# example scenario
ls_df = scale_load_shed_results(224, 1998)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls)

if_flow_df = scale_if_flow_results(224, 1998)

ax1 = Axis(fig[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=-0.5:0.5:1.0,
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
           title="(c) NYS-PJM interface utilization",
           titlegap=5.0
)
vspan!(ax3, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax3, if_flow_df[15, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

linkxaxes!(ax1, ax3)
linkxaxes!(ax2, ax3)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace
ax3.yticklabelspace = yspace

colsize!(fig.layout, 1, Aspect(1, 12.0))
rowgap!(fig.layout, 10.0)
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "clus1_ABif_BCif_PJMif_utilization.png"), fig)