using Distributed

num_procs = parse(Int64, ENV["SLURM_NTASKS"])
addprocs(num_procs, topology=:master_worker)

# activate environment
@everywhere begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", ".."))
    Pkg.instantiate()
end

# load dependencies and utility functions
@everywhere begin
    using CairoMakie
    using CSV
    using DataFrames
    using DelimitedFiles
    using ProgressMeter
    using Random

    include(joinpath(@__DIR__, "..", "utils", "kmodes.jl"))
    include(joinpath(@__DIR__, "..", "utils", "other.jl"))

    Random.seed!(myid()) # for reproducibility

    data_dir = joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "independence")
    output_dir = joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "clustering")
end

# main script
# -----------

# collect results from independence tests
df = DataFrame(vcat(:scenario => repeat(1:300, outer=22), 
               :yr => repeat(1998:2019, inner=300), 
               [Symbol("if$(i)") => zeros(Bool, 6600) for i in 1:15])
)
for yr in 1998:2019
    tmp_df = DataFrame(CSV.File(joinpath(data_dir, "if_flow_res_$(yr).csv")))
    gdf = groupby(tmp_df, :scenario)
    df[df.yr .== yr, Between(:if1, :if15)] = mapreduce(x -> (x.bus_44_pvalue .< 0.05)', vcat, gdf) # binary encoding
end

# 0/1 matrices where 1 = significant nonzero TE, 0 = no significant nonzero TE
data = Matrix(df[:, Not(:scenario, :yr)])'

# try k-modes clustering with 2 ≤ k ≤ 20
cluster_res_df = DataFrame(k=2:20, 
                           obj=Vector{Float64}(undef, 19), 
                           n_iter=Vector{Int}(undef, 19), 
                           min_cluster_size=Vector{Int}(undef, 19), 
                           max_cluster_size=Vector{Int}(undef, 19)
)
n_init = 300 # try different initializations
for (j, k) in enumerate(2:20)
    res = @showprogress pmap(1:n_init) do i
        kmodes(data, k)
    end

    opt_res = res[argmin([r.total_cost for r in res])]
    cluster_res_df[j, :obj] = opt_res.total_cost
    cluster_res_df[j, :n_iter] = opt_res.iterations
    cluster_res_df[j, :min_cluster_size] = minimum(opt_res.counts)
    cluster_res_df[j, :max_cluster_size] = maximum(opt_res.counts)

    writedlm(joinpath(output_dir, "cluster_labels_$(k).txt"), opt_res.assignments)
    writedlm(joinpath(output_dir, "cluster_modes_$(k).txt"), opt_res.modes)
end

# save summary results
CSV.write(joinpath(output_dir, "cluster_res.csv"), cluster_res_df)