using Distributed

num_procs = parse(Int64, ENV["SLURM_NTASKS"])
addprocs(num_procs, topology=:master_worker)

# activate environment
@everywhere begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", ".."))
    Pkg.instantiate()
end

# load dependencies and utility functions
@everywhere begin
    using ProgressMeter
    using Random

    include(joinpath(@__DIR__, "setup.jl"))

    Random.seed!(myid()) # for reproducibility
    global const JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))
end

# main script
# -----------

yrs = 1998:2019
yr = yrs[parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])]
load_bus = 44 # bus 44 is only target

# independence test configuration
cfg = IndependenceTestConfig(k_embed=50, k_est=20, bins=range(0.0, 1.2, length=7), nshuffles=200, w=23)

for scenario in 1:300
    # to store results
    if_flow_res_df = DataFrame("id" => 1:15,
                               "scenario" => scenario*ones(Int, 15),
                               "bus_$(load_bus)_pvalue" => zeros(Float64, 15),
                               "bus_$(load_bus)_est" => zeros(Float64, 15),
                               "bus_$(load_bus)_min" => zeros(Float64, 15),
                               "bus_$(load_bus)_max" => zeros(Float64, 15)
    )

    # load target data
    ls_df = scale_load_shed_results(scenario, yr)
    s_trg = ls_df[load_bus, :scaled_load_shed][JJA_idx]

    # optimize target embedding
    opt = EmbeddingOptions(max_d=8, max_τ=24, nbhd_type="mass", nbhd_param=cfg.k_embed, w=23, metric=Cityblock())
    dT, τT = optimize_embedding(s_trg, opt; n_itr=20)
    println("Scenario $(scenario): dT = $(dT), τT = $(τT)")

    # perform independence tests
    res = @showprogress pmap(1:nrow(if_flow_res_df)) do i
        if_flow_independence(s_trg, i, scenario, yr, cfg, dT, τT)
    end

    # save results
    if_flow_res_df[:, Symbol("bus_$(load_bus)_pvalue")] = [r.pvalue for r in res]
    if_flow_res_df[:, Symbol("bus_$(load_bus)_est")]    = [r.m for r in res]
    if_flow_res_df[:, Symbol("bus_$(load_bus)_min")]    = [minimum(r.m_surr) for r in res]
    if_flow_res_df[:, Symbol("bus_$(load_bus)_max")]    = [maximum(r.m_surr) for r in res]
    
    outfile = joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "independence", "if_flow_res_$(yr).csv")
    if isfile(outfile)
        CSV.write(outfile, if_flow_res_df, append=true)
    else
        CSV.write(outfile, if_flow_res_df)
    end
end