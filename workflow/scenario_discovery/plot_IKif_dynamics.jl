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
Jl_idx = in(get_month_idx(7)).(JJA_idx)

# scenario_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "output", "scenario_features.csv")))
# labels = readdlm(joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "cluster_labels_3.txt"), Int)
# insertcols!(scenario_df, :label => vec(labels))

# clus1_ids = vec(labels .== 1)
# clus1_df = DataFrame(:scenario => scenario_df[clus1_ids, :scenario],
#                      :yr => scenario_df[clus1_ids, :yr],
#                      :prop_IK_congestion_zoneGls1 => Vector{Float64}(undef, sum(clus1_ids)), # I-K congestion with zone G load shedding
#                      :prop_IK_congestion_zoneGls0 => Vector{Float64}(undef, sum(clus1_ids)), # I-K congestion without zone G load shedding
#                      :avg_if9_flow => Vector{Float64}(undef, sum(clus1_ids)), # average H-I flow
#                      :avg_if9_flow_zoneJKls1 => Vector{Float64}(undef, sum(clus1_ids)), # average H-I flow with zone J/K load shedding
#                      :avg_if9_util => Vector{Float64}(undef, sum(clus1_ids)), # average H-I utilization
#                      :avg_if9_util_zoneJKls1 => Vector{Float64}(undef, sum(clus1_ids)) # average H-I utilization with zone J/K load shedding
# )
# for i in 1:nrow(clus1_df)
#     scenario, yr = clus1_df.scenario[i], clus1_df.yr[i]
    
#     ls_df = scale_load_shed_results(scenario, yr)
#     lsG = ls_df[44, :load_shed][JJA_idx]
#     lsJK = sum(ls_df[46:49, :load_shed])[JJA_idx]

#     if_flow_df = scale_if_flow_results(scenario, yr)
#     if11_util = if_flow_df[11, :scaled_if_flow][JJA_idx]
#     clus1_df[i, :prop_IK_congestion_zoneGls1] = sum(if11_util[lsG .> 0.01] .>= 0.9) / sum(abs.(if11_util[lsG .> 0.01]) .>= 0.9)
#     clus1_df[i, :prop_IK_congestion_zoneGls0] = sum(if11_util[lsG .<= 0.01] .>= 0.9) / sum(abs.(if11_util[lsG .<= 0.01]) .>= 0.9)

#     if9_flow = if_flow_df[9, :if_flow][JJA_idx]
#     clus1_df[i, :avg_if9_flow] = mean(if9_flow)
#     clus1_df[i, :avg_if9_flow_zoneJKls1] = mean(if9_flow[lsJK .> 0.01])

#     if9_util = if_flow_df[9, :scaled_if_flow][JJA_idx]
#     clus1_df[i, :avg_if9_util] = mean(if9_util)
#     clus1_df[i, :avg_if9_util_zoneJKls1] = mean(if9_util[lsJK .> 0.01])
# end
# CSV.write(joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "clus1_df.csv"), clus1_df)

# read in data (or run above block)
clus1_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "clus1_df.csv")))

# plot 
pt = 4/3
inch = 96
fig = Figure(size=(0inch, 3.0inch), fontsize=10pt)
hex_colors = [
    "#a7e4d3",  # Pale turquoise
    "#6fd3b7",  # Light sea green
    "#38c29b",  # Mid-tone teal
    "#019e73",  # Base green
    "#017b5a",  # Deep teal
    "#015640",  # Very dark green
    "#003b2c"   # Almost black green
]
cmap = parse.(Colorant, hex_colors)
p_if9 = sortperm(clus1_df.avg_if9_util, rev=true)

gright = GridLayout()
axright = gright[1,1] = Axis(fig,
                             xgridvisible=false, ygridvisible=false,
                             xlabel="without zone G power shortage", ylabel="with zone G power shortage", 
                             limits=((0.0, 0.5), (0.0, 1.0)),
                             xticks=[0.0, 0.5],
                             title="(c) Prop. of elevated flow I → K", titlegap=8.0
)
poly!(axright, Point2f[(0,0.5), (0.5,0.5), (0.5,1.0), (0,1.0)], color="#f0f0f0")
sc = scatter!(axright, 
              clus1_df.prop_IK_congestion_zoneGls0[p_if9], clus1_df.prop_IK_congestion_zoneGls1[p_if9], 
              color=clus1_df.avg_if9_util[p_if9], colormap=cmap, alpha=0.8
)
ids = [findall((clus1_df.scenario .== 24) .&& (clus1_df.yr .== 2011))[1]; findall((clus1_df.scenario .== 115) .&& (clus1_df.yr .== 2005))]
scatter!(axright, 
         clus1_df.prop_IK_congestion_zoneGls0[ids], clus1_df.prop_IK_congestion_zoneGls1[ids],
         color=clus1_df.avg_if9_util[ids], colormap=cmap, strokecolor=:black, strokewidth=1.5
)
text!(axright, 
      clus1_df.prop_IK_congestion_zoneGls0[ids], clus1_df.prop_IK_congestion_zoneGls1[ids],
      text=["A", "B"], align=(:left, :bottom), offset=(2.0, 2.0)
)
Colorbar(fig[1,3], 
         sc, 
         label="mean H-I interface utilization", 
         vertical=true
)

gleft = GridLayout()
axtop = gleft[1,1] = Axis(fig,
                          xgridvisible=false, ygridvisible=false,
                          xlabel="hour", ylabel="utilization",
                          limits=((1, sum(Jl_idx)), nothing),
                          yticks=-1.0:1.0:1.0, 
                          ytickformat="{:.1f}",
                          title="(a) I-K interface utilization in Scenario A", titlegap=5.0
)
hidexdecorations!(axtop)
ls_df = scale_load_shed_results(24, 2011)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls[Jl_idx])
vspan!(axtop, xlow, xhigh, color=("#d32f2f", 0.3))
if_flow_df = scale_if_flow_results(24, 2011)
lines!(axtop, if_flow_df[11, :scaled_if_flow][JJA_idx][Jl_idx], color=:gray22, linewidth=1.0)

axbottom = gleft[2,1] = Axis(fig,
                             xgridvisible=false, ygridvisible=false,
                             xlabel="hour", ylabel="utilization",
                             limits=((1, sum(Jl_idx)), nothing),
                             xticks=vcat(1, 168:168:sum(Jl_idx)), yticks=-1.0:1.0:1.0,
                             ytickformat="{:.1f}",
                             title="(b) I-K interface utilization in Scenario B", titlegap=5.0
)
ls_df = scale_load_shed_results(115, 2005)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls[Jl_idx])
vspan!(axbottom, xlow, xhigh, color=("#d32f2f", 0.3))
if_flow_df = scale_if_flow_results(115, 2005)
lines!(axbottom, if_flow_df[11, :scaled_if_flow][JJA_idx][Jl_idx], color=:gray22, linewidth=1.0)

fig.layout[1,1] = gleft
fig.layout[1,2] = gright
colsize!(fig.layout, 2, Aspect(1, 0.5))
colsize!(fig.layout, 1, Aspect(1, 2.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "clus1_IKif_util.png"), fig)
