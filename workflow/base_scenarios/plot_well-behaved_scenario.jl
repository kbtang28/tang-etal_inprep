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
Ju_idx = in(get_month_idx(6)).(JJA_idx)
JJ_idx = in(vcat(get_month_idx(6), get_month_idx(7))).(JJA_idx)

pt = 4/3
inch = 96

# plot G-H and F-G IF utilization, zone F renewables availability for well-behaved scenario
fig = Figure(size=(0inch, 2.25inch), fontsize=10pt)
_, _, _, wind_cap_sf, solar_cap_sf, _, du_scenario = get_du_factors(140)
ls_df = scale_load_shed_results(140, 2012)
ls = ls_df[44, :scaled_load_shed][JJA_idx]
xlow, xhigh = get_ls_span(ls[Jl_idx])

if_flow_df = scale_if_flow_results(140, 2012)

ax1 = Axis(fig[1,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,sum(Jl_idx)), nothing),
           xticks=vcat(1, 168:168:sum(Jl_idx)),
           xlabel="hour", ylabel="utilization",
           title="(a) G-H interface utilization",
           titlegap=6.0
)
hidexdecorations!(ax1)
vspan!(ax1, xlow, xhigh, color=("#d32f2f", 0.3), label="hrs with power shortage")
lines!(ax1, if_flow_df[8, :scaled_if_flow][JJA_idx][Jl_idx], color=:gray22, linewidth=1.0)

wind_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "wind", "wind_2012.csv"), ',', Float64)[:, 2:end]
wind_cap_mat = wind_cap_sf*round.(wind_cap_mat[:, JJA_idx]; digits=2)
wind_cap_F = wind_cap_mat[12, :]
solar_cap_mat = readdlm(joinpath(@__DIR__, "..", "..", "ACORN", "solar", "Scenario$(du_scenario)", "solarUPV2012.csv"), ',', Float64)[:, 2:end]
solar_cap_mat = solar_cap_sf*round.(solar_cap_mat[:, JJA_idx]; digits=2)
solar_cap_F = solar_cap_mat[10, :]
ax2 = Axis(fig[2,1],
           xgridvisible=false, ygridvisible=false,
           limits=((1,sum(Jl_idx)), nothing),
           xticks=vcat(1, 168:168:sum(Jl_idx)),
           xlabel="hour", ylabel="availability\n(10⁴ MWh)",
           ytickformat="{:.1f}",
           title="(b) Zone F renewable energy availability",
           titlegap=6.0
)
hidexdecorations!(ax2)
vspan!(ax2, xlow, xhigh, color=("#d32f2f", 0.3))
band!(ax2, 
      1:sum(Jl_idx), zeros(sum(Jl_idx)), wind_cap_F[Jl_idx]/1e4, 
      color="#99d8c7", label="wind"
)
band!(ax2, 
      1:sum(Jl_idx), wind_cap_F[Jl_idx]/1e4, (wind_cap_F[Jl_idx].+solar_cap_F[Jl_idx])/1e4, 
      color="#eec680", label="solar"
)
lines!(ax2, (wind_cap_F[Jl_idx].+solar_cap_F[Jl_idx])/1e4, color=:gray22, linewidth=1.0)
axislegend(ax2, patchsize=(10.0,10.0), position=:lt, labelsize=8pt, padding=(3.0,3.0,3.0,3.0))

linkxaxes!(ax1, ax2)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace+3.0
ax2.yticklabelspace = yspace+3.0

colsize!(fig.layout, 1, Aspect(1, 8.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "..", "figures", "s140_GHif_zoneFrenewables.png"), fig)