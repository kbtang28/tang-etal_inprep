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
Jl_idx = in(get_month_idx(7)).(JJA_idx)

pt = 4/3
inch = 96

# plot A-B, C-E, G-H, H-I, I-K IF utilization for extreme temperature scenario
fig = Figure(size=(0inch, 5.0inch), fontsize=10pt)
ls_df = scale_load_shed_results(290, 2002)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls)

if_flow_df = scale_if_flow_results(290, 2002)

ax1 = Axis(fig[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.5:1.0,
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
           title="(b) C-E interface utilization",
           titlegap=5.0
)
hidexdecorations!(ax2)
vspan!(ax2, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax2, if_flow_df[3, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

ax3 = Axis(fig[3,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)),
           xlabel="hour", ylabel="utilization",
           title="(c) G-H interface utilization",
           titlegap=5.0
)
hidexdecorations!(ax3)
vspan!(ax3, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax3, if_flow_df[8, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

ax4 = Axis(fig[4,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=0.0:0.4:0.8,
           xlabel="hour", ylabel="utilization",
           title="(d) H-I interface utilization",
           titlegap=5.0
)
hidexdecorations!(ax4)
vspan!(ax4, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax4, if_flow_df[9, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

ax5 = Axis(fig[5,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,length(JJA_idx)), nothing),
           xticks=vcat(1, 168:168:length(JJA_idx)), yticks=-1.0:1.0:1.0,
           ytickformat="{:.1f}",
           xlabel="hour", ylabel="utilization",
           title="(e) I-K interface utilization",
           titlegap=5.0
)
vspan!(ax5, xlow, xhigh, color=("#d32f2f", 0.3))
lines!(ax5, if_flow_df[11, :scaled_if_flow][JJA_idx], color=:gray22, linewidth=1.0)

linkxaxes!(ax1, ax5)
linkxaxes!(ax2, ax5)
linkxaxes!(ax3, ax5)
linkxaxes!(ax4, ax5)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3, ax4, ax5])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace
ax3.yticklabelspace = yspace
ax4.yticklabelspace = yspace
ax5.yticklabelspace = yspace

colsize!(fig.layout, 1, Aspect(1, 12.0))
rowgap!(fig.layout, 10.0)
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "s290_ABif_CEif_GHif_HIif_IKif_utilization.png"), fig)

# plot H-I, I-K IF utilization (with scatterplot) for extreme temperature scenario
fig = Figure(size=(8.0inch, 2.75inch), fontsize=10pt)
ls_df = scale_load_shed_results(290, 2002)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls[Jl_idx])

if_flow_df = scale_if_flow_results(290, 2002)
if9 = if_flow_df[9, :scaled_if_flow][JJA_idx]
if11 = if_flow_df[11, :scaled_if_flow][JJA_idx]

gleft = GridLayout()
axleft = gleft[1,1] = Axis(fig,
                           xgridvisible=false, ygridvisible=false,
                           xlabel="Interface H-I utilization",
                           ylabel="Interface I-K utilization"
)
hlines!(axleft, [0.0], color=:black, linestyle=:dash)
scatter!(axleft, if9, if11, color=("#7f7f7f", 0.3))
text!(axleft, [(-0.05, 0.5)], color="#7f7f7f", align=(:center, :center), rotation=pi/2, text="I → K flow")
text!(axleft, [(-0.05, -0.5)], color="#7f7f7f", align=(:center, :center), rotation=pi/2, text="K → I flow")
text!(axleft, [(maximum(if9), minimum(if11))], text="(c)", font=:bold, align=(:right, :bottom))
gright = GridLayout()
axtop = gright[1,1] = Axis(fig,
                           xgridvisible=false, ygridvisible=false,
                           xlabel="hour", ylabel="utilization",
                           limits=((1, sum(Jl_idx)), nothing),
                           yticks=-1.0:1.0:1.0,
                           title="(a) I-K interface utilization"
)
hidexdecorations!(axtop)
vspan!(axtop, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(axtop, if11[Jl_idx], color=:gray22, linewidth=1.0)

axbottom = gright[2,1] = Axis(fig,
                              xgridvisible=false, ygridvisible=false,
                              xlabel="hour", ylabel="utilization",
                              limits=((1, sum(Jl_idx)), nothing),
                              xticks=vcat(1, 168:168:sum(Jl_idx)),
                              title="(b) H-I interface utilization"
)
vspan!(axbottom, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(axbottom, if9[Jl_idx], color=:gray22, linewidth=1.0)

fig.layout[1,2] = gleft
fig.layout[1,1] = gright
colsize!(fig.layout, 2, Aspect(1, 1.0))
rowgap!(gright, 1, 5.0)
save(joinpath(@__DIR__, "..", "..", "figures", "s290_HIif_IKif_utilization.png"), fig)
