# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# load dependencies and helper functions
using CairoMakie
using DelimitedFiles

import ColorSchemes: seaborn_colorblind

include(joinpath(@__DIR__, "utils", "scaling.jl"))

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))
Jl_idx = 721:1464 # subset of JJA_idx representing month of July

# plot
pt = 4/3
inch = 96
fig = Figure(size=(0inch, 2.75inch), fontsize=10pt)

# interface example
if_df = scale_if_flow_results(140, 2012)
s = if_df[12, :scaled_if_flow][JJA_idx]
surros = readdlm(joinpath(@__DIR__, "base_scenarios", "precomputed_wls_surrogates", "if_flow12_surro_s140_2012.txt"))
ax1 = Axis(fig[1,1], 
           xgridvisible=false, ygridvisible=false, 
           limits=((1, length(Jl_idx)), nothing), 
           title="Interface example", 
           xlabel="hour", ylabel="utilization", 
           xticks=(vcat(1, 168:168:length(Jl_idx)), string.(vcat(1, 168:168:length(Jl_idx))))
)
hidexdecorations!(ax1)
for i in 1:10
    lines!(ax1, 2*surros[Jl_idx, i], color=(:gray69, 0.4), linewidth=0.75)
end
lines!(ax1, s[Jl_idx], color=seaborn_colorblind[1])

# curtailment example
curtail_df = scale_curtail_results(140, 2012)
s = curtail_df[26, :scaled_curtail][JJA_idx]
surros = readdlm(joinpath(@__DIR__, "base_scenarios", "precomputed_wls_surrogates", "curtail26_surro_s140_2012.txt"))
ax2 = Axis(fig[2,1], 
           xgridvisible=false, ygridvisible=false, 
           limits=((1, length(Jl_idx)), nothing), 
           title="Solar generator example", 
           xlabel="hour", ylabel="curtailment", 
           xticks=(vcat(1, 168:168:length(Jl_idx)), string.(vcat(1, 168:168:length(Jl_idx))))
)
for i in 1:10
    lines!(ax2, surros[Jl_idx, i], color=(:gray69, 0.4), linewidth=0.75)
end
lines!(ax2, s[Jl_idx], color=seaborn_colorblind[1])

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace+5.0
ax2.yticklabelspace = yspace+5.0

colsize!(fig.layout, 1, Aspect(1, 8.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "figures", "wls_surro_examples.png"), fig, dpi=300)