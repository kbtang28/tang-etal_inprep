# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# load dependencies and helper functions
using CairoMakie
using CSV
using DataFrames

import ColorSchemes: seaborn_colorblind

include(joinpath(@__DIR__, "utils", "scaling.jl"))

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))

# ls_hrs_df = DataFrame(Dict(Symbol("bus$(i)") => Vector{Float64}(undef, 6600) for i in 1:57))
# ls_hrs_df.scenario = repeat(1:300, outer=22)
# ls_hrs_df.yr = repeat(1998:2019, inner=300)

# ls_prop_df = DataFrame(Dict(Symbol("bus$(i)") => 100.0*ones(Float64, 6600) for i in 1:57))
# ls_prop_df.scenario = repeat(1:300, outer=22)
# ls_prop_df.yr = repeat(1998:2019, inner=300)

# for i in 1:6600
#     scenario, yr = ls_hrs_df[i, [:scenario, :yr]]
#     ls_df = scale_load_shed_results(scenario, yr)
#     total_load = sum(calc_load(scenario, yr)[:, JJA_idx], dims=2)
    
#     for b in 1:57
#         ls = ls_df[b, :load_shed][JJA_idx]
#         ls_hrs_df[i, Symbol("bus$(b)")] = sum(ls .> 0.0)
#         if total_load[b] > 0.0
#             ls_prop_df[i, Symbol("bus$(b)")] = sum(ls)/total_load[b]
#         end
#     end
# end

# CSV.write(joinpath(@__DIR__, "..", "output", "ls_hrs_by_bus.csv"), ls_hrs_df)
# CSV.write(joinpath(@__DIR__, "..", "output", "ls_prop_by_bus.csv"), ls_prop_df)

# read in data (or run above block)
ls_hrs_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "output", "ls_hrs_by_bus.csv")))
ls_prop_df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "output", "ls_prop_by_bus.csv")))

zone_dict = Dict("A" => [21, 22, 23, 24, 25, 26, 27, 28],
                 "B" => [19, 20, 29],
                 "C" => [17, 18, 30, 31, 32, 33, 34, 35, 37, 38, 39],
                 "D" => [15, 16],
                 "E" => [5, 10, 11, 12, 13, 14, 36],
                 "F" => [4, 7, 8, 9],
                 "G" => [6, 40, 42, 43, 44],
                 "H" => [41],
                 "I" => [45],
                 "J" => [48, 49],
                 "K" => [46, 47],
                 "PJM" => [53, 54, 55, 56, 57],
                 "NE" => [1, 2, 3],
                 "IESO" => [50, 51, 52]
)
NYS_ids = 4:49 # buses 1-3 belong to ISONE, buses 50-52 to IESO, buses 53-57 to PJM
bus_dict = Dict(bus => zone for (zone, buses) in zone_dict for bus in buses)

# plot
pt = 4/3
inch = 96
fig = Figure(size=(3.5inch, 3.5inch), fontsize=10pt)

ids = intersect(findall([sum(ls_hrs_df[:, Symbol("bus$(b)")]) for b in 1:57] .> 0.0), NYS_ids)
zones = [bus_dict[id] for id in ids]
ids = ids[sortperm(zones)]
id44 = findfirst(44 .== ids)

new_ls_hrs_df = stack(ls_hrs_df[:, [Symbol("bus$(b)") for b in ids]], variable_name=:bus_id, value_name=:hrs)
transform!(new_ls_hrs_df, :bus_id => ByRow(str -> parse(Int64, strip(str, ['b', 'u', 's']))) => :bus_id)
transform!(new_ls_hrs_df, :bus_id => ByRow(id -> findfirst(id .== ids)) => :bus_id)

ax1 = Axis(fig[1,1],
           xgridvisible=false, ygridvisible=false,
           xlabel="bus zone", ylabel="prop. of hours\nwith unmet demand",
           xticks=(1:length(ids), zones[sortperm(zones)])
)
hidexdecorations!(ax1)
boxplot!(ax1, 
         new_ls_hrs_df[new_ls_hrs_df.bus_id .!= id44, :bus_id], new_ls_hrs_df[new_ls_hrs_df.bus_id .!= id44, :hrs]/length(JJA_idx),
         color=seaborn_colorblind[1]
)
boxplot!(ax1,
         new_ls_hrs_df[new_ls_hrs_df.bus_id .== id44, :bus_id], new_ls_hrs_df[new_ls_hrs_df.bus_id .== id44, :hrs]/length(JJA_idx),
         color=seaborn_colorblind[2]
)
text!(11, maximum(new_ls_hrs_df.hrs)/length(JJA_idx), text="(a)", font=:bold, align=(:left, :top))

new_ls_prop_df = stack(ls_prop_df[:, [Symbol("bus$(b)") for b in ids]], variable_name=:bus_id, value_name=:ls_prop)
transform!(new_ls_prop_df, :bus_id => ByRow(str -> parse(Int64, strip(str, ['b', 'u', 's']))) => :bus_id)
transform!(new_ls_prop_df, :bus_id => ByRow(id -> findfirst(id .== ids)) => :bus_id)

ax2 = Axis(fig[2,1],
           xgridvisible=false, ygridvisible=false,
           xlabel="bus zone", ylabel="prop. of demand unmet",
           xticks=(1:length(ids), zones[sortperm(zones)])
)
boxplot!(ax2, 
         new_ls_prop_df[new_ls_prop_df.bus_id .!= id44, :bus_id], new_ls_prop_df[new_ls_prop_df.bus_id .!= id44, :ls_prop],
         color=seaborn_colorblind[1]
)
boxplot!(ax2,
         new_ls_prop_df[new_ls_prop_df.bus_id .== id44, :bus_id], new_ls_prop_df[new_ls_prop_df.bus_id .== id44, :ls_prop],
         color=seaborn_colorblind[2]
)
text!(11, maximum(new_ls_prop_df.ls_prop), text="(b)", font=:bold, align=(:left, :top))

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace+2.0
ax2.yticklabelspace = yspace+3.0

colsize!(fig.layout, 1, Aspect(1, 2.0))
resize_to_layout!(fig)

save(joinpath(@__DIR__, "..", "figures", "ls_hrs_and_prop_by_bus.png"), fig)