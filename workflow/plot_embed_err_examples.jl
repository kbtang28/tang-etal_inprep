# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# load dependencies and helper functions
using CairoMakie
using CSV
using DataFrames

import ColorSchemes: amp

include(joinpath(@__DIR__, "utils", "scaling.jl"))
include(joinpath(@__DIR__, "utils", "embedding.jl"))

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))
Jl_idx = 721:1464 # month of July

# load data
if_flow_df = scale_if_flow_results(140, 2012)

# more periodic example
s11 = if_flow_df[11, :scaled_if_flow][JJA_idx]
opt = EmbeddingOptions(max_d=8, max_τ=24, nbhd_type="mass", nbhd_param=50, w=23, metric=Cityblock())
query_idx = collect(1+(opt.max_d-1)*opt.max_τ+opt.u : length(s11))
mre11 = zeros(opt.max_τ, opt.max_d)
for τ in 1:opt.max_τ, d in 1:opt.max_d
    for ct in 1:10
        local_pred = find_local_pred(s11, query_idx, d, τ, opt)
        mre11[τ, d] += ragwitz_mre(s11, query_idx, local_pred)
    end
end
mre11 ./= 10

# less periodic example
s14 = if_flow_df[14, :scaled_if_flow][JJA_idx]
mre14 = zeros(opt.max_τ, opt.max_d)
for τ in 1:opt.max_τ, d in 1:opt.max_d
    for ct in 1:10
        local_pred = find_local_pred(s14, query_idx, d, τ, opt)
        mre14[τ, d] += ragwitz_mre(s14, query_idx, local_pred)
    end
end
mre14 ./= 10

# plot
pt = 4/3
inch = 96
fig = Figure(size=(7inch, 3.5inch), fontsize=10pt)

ax11 = Axis(fig[1,1],
            limits=((1,336), nothing),
            ylabel="utilization",
            xlabel="hour",
            xticks=vcat(1, 48:48:336),
            xgridvisible=false,
            ygridvisible=false,
            title="NY-IESO interface",
)
ax21 = Axis(fig[2,1],
            limits=((1,336), nothing),
            ylabel="utilization",
            xlabel="hour",
            xticks=vcat(1, 48:48:336),
            xgridvisible=false,
            ygridvisible=false,
            title="I-K interface",
)
ax12 = Axis(fig[1,2], 
            ylabel="dimension", 
            xlabel="delay", 
            yticks=2:2:opt.max_d, 
            xticks=2:2:opt.max_τ, 
            title="MRSE for local model predictions", 
)
ax22 = Axis(fig[2,2], 
            ylabel="dimension", 
            xlabel="delay", 
            yticks=2:2:opt.max_d, 
            xticks=2:2:opt.max_τ, 
            title="MRSE for local model predictions", 
)
heatmap!(ax12, 1:opt.max_τ, 1:opt.max_d, mre14, colormap=:amp)
scatter!(ax12, findmin(mre14)[2][1], findmin(mre14)[2][2], marker=:xcross, color=:black, markersize=:10)

heatmap!(ax22, 1:opt.max_τ, 1:opt.max_d, mre11, colormap=:amp)
scatter!(ax22, findmin(mre11)[2][1], findmin(mre11)[2][2], marker=:xcross, color=:black, markersize=:10)

Colorbar(fig[1, 3], colormap=:amp, colorrange=extrema(mre14), tickformat="{:.2f}")
Colorbar(fig[2, 3], colormap=:amp, colorrange=extrema(mre11), tickformat="{:.2f}")

lines!(ax11, s14[Jl_idx[1:336]], linewidth=1.0, color=:grey22)
s14_est = vcat(zeros(length(JJA_idx)-length(query_idx)), find_local_pred(s14, query_idx, 4, 11, opt)) # predictions with optimal embedding parameters
lines!(ax11, s14_est[Jl_idx[1:336]], linewidth=1.0, color=("#d32f2f", 0.5))
lines!(ax21, s11[Jl_idx[1:336]], linewidth=1.0, color=:grey22)
s11_est = vcat(zeros(length(JJA_idx)-length(query_idx)), find_local_pred(s11, query_idx, 3, 1, opt)) # predictions with optimal embedding parameters
lines!(ax21, s11_est[Jl_idx[1:336]], linewidth=1.0, color=("#d32f2f", 0.5))

colsize!(fig.layout, 1, Aspect(1, 2.5))
colsize!(fig.layout, 2, Aspect(1, 3.0))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "..", "figures", "ragwitz_criterion.png"), fig)