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
    using Associations
    using DataFrames
    using ProgressMeter

    include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))
    include(joinpath(@__DIR__, "..", "utils", "independence_tests.jl"))
    include(joinpath(@__DIR__, "setup.jl"))

    Random.seed!(myid()) # for reproducibility
    global const JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))
end

# main script
# -----------

wavelet_surros_dir = joinpath(@__DIR__, "precomputed_wls_surrogates")
output_dir = joinpath(@__DIR__, "..", "..", "output", "base_scenarios")

scenario_yr_pairs = [(140, 2012), (290, 2002), (69, 2011)]
scenario, yr = scenario_yr_pairs[parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])]
load_bus = 44

# to store results
batt_df = scale_batt_results(scenario, yr)
select!(batt_df, Not(:state, :scaled_state))

# to store results
branches = [99, 100] # last two branches are new HVDC lines
flow_df = scale_flow_results(scenario, yr)[branches, :]
select!(flow_df, Not(:flow, :scaled_flow))

# to store results
if_flow_df = scale_if_flow_results(scenario, yr)
select!(if_flow_df, Not(:if_flow, :scaled_if_flow))

# to store results
gens = 1:29 # dispatchable energy sources (hydro & imported)
gen_df = scale_gen_results(scenario, yr)[gens, :]
select!(gen_df, Not(:gen, :scaled_gen))

# to store results
curtail_df = scale_curtail_results(scenario, yr)
select!(curtail_df, Not(:curtail, :scaled_curtail))

# test batteries - wavelet surrogates
res = @showprogress pmap(1:nrow(batt_df)) do i
    surro_file = joinpath(wavelet_surros_dir, "batt$(i)_surro_s$(scenario)_$(yr).txt")
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)
    
    return batt_independence(load_bus, i, scenario, yr, cfg, PrecomputedWLS(file=surro_file))
end

# save results
insertcols!(batt_df, Symbol("wls_pvalue") => [r.pvalue for r in res])
insertcols!(batt_df, Symbol("wls_est") => [r.m for r in res])
insertcols!(batt_df, Symbol("wls_min") => [minimum(r.m_surr) for r in res])
insertcols!(batt_df, Symbol("wls_max") => [maximum(r.m_surr) for r in res])

# test batteries - block-shuffled surrogates
res = @showprogress pmap(1:nrow(batt_df)) do i
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return batt_independence(load_bus, i, scenario, yr, cfg, BlockShuffle(2; shift=true))
end

# save results
insertcols!(batt_df, Symbol("block_pvalue") => [r.pvalue for r in res])
insertcols!(batt_df, Symbol("block_est") => [r.m for r in res])
insertcols!(batt_df, Symbol("block_min") => [minimum(r.m_surr) for r in res])
insertcols!(batt_df, Symbol("block_max") => [maximum(r.m_surr) for r in res])

# test flows - wavelet surrogates
res = @showprogress pmap(flow_df.id) do i
    surro_file = joinpath(wavelet_surros_dir, "flow$(i)_surro_s$(scenario)_$(yr).txt")
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return flow_independence(load_bus, i, scenario, yr, cfg, PrecomputedWLS(file=surro_file))
end

# save results
insertcols!(flow_df, Symbol("wls_pvalue") => [r.pvalue for r in res])
insertcols!(flow_df, Symbol("wls_est") => [r.m for r in res])
insertcols!(flow_df, Symbol("wls_min") => [minimum(r.m_surr) for r in res])
insertcols!(flow_df, Symbol("wls_max") => [maximum(r.m_surr) for r in res])

# test flows - block-shuffled surrogates
res = @showprogress pmap(flow_df.id) do i
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return flow_independence(load_bus, i, scenario, yr, cfg, BlockShuffle(2; shift=true))
end

# save results
insertcols!(flow_df, Symbol("block_pvalue") => [r.pvalue for r in res])
insertcols!(flow_df, Symbol("block_est") => [r.m for r in res])
insertcols!(flow_df, Symbol("block_min") => [minimum(r.m_surr) for r in res])
insertcols!(flow_df, Symbol("block_max") => [maximum(r.m_surr) for r in res])

# test IF flows - wavelet surrogates
res = @showprogress pmap(1:nrow(if_flow_df)) do i
    surro_file = joinpath(wavelet_surros_dir, "if_flow$(i)_surro_s$(scenario)_$(yr).txt")
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return if_flow_independence(load_bus, i, scenario, yr, cfg, PrecomputedWLS(file=surro_file))
end

# save results
insertcols!(if_flow_df, Symbol("wls_pvalue") => [r.pvalue for r in res])
insertcols!(if_flow_df, Symbol("wls_est") => [r.m for r in res])
insertcols!(if_flow_df, Symbol("wls_min") => [minimum(r.m_surr) for r in res])
insertcols!(if_flow_df, Symbol("wls_max") => [maximum(r.m_surr) for r in res])

# test IF flows - block-shuffled surrogates
res = @showprogress pmap(1:nrow(if_flow_df)) do i
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return if_flow_independence(load_bus, i, scenario, yr, cfg, BlockShuffle(2; shift=true))
end

# save results
insertcols!(if_flow_df, Symbol("block_pvalue") => [r.pvalue for r in res])
insertcols!(if_flow_df, Symbol("block_est") => [r.m for r in res])
insertcols!(if_flow_df, Symbol("block_min") => [minimum(r.m_surr) for r in res])
insertcols!(if_flow_df, Symbol("block_max") => [maximum(r.m_surr) for r in res])

# test generators - wavelet surrogates
res = @showprogress pmap(1:nrow(gen_df)) do i
    surro_file = joinpath(wavelet_surros_dir, "gen$(i)_surro_s$(scenario)_$(yr).txt")
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)
    
    return gen_independence(load_bus, i, scenario, yr, cfg, PrecomputedWLS(file=surro_file))
end

# save results
insertcols!(gen_df, Symbol("wls_pvalue") => [r.pvalue for r in res])
insertcols!(gen_df, Symbol("wls_est") => [r.m for r in res])
insertcols!(gen_df, Symbol("wls_min") => [minimum(r.m_surr) for r in res])
insertcols!(gen_df, Symbol("wls_max") => [maximum(r.m_surr) for r in res])

# test generators - block-shuffled surrogates
res = @showprogress pmap(1:nrow(gen_df)) do i
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)
    
    return gen_independence(load_bus, i, scenario, yr, cfg, BlockShuffle(2; shift=true))
end

# save results
insertcols!(gen_df, Symbol("block_pvalue") => [r.pvalue for r in res])
insertcols!(gen_df, Symbol("block_est") => [r.m for r in res])
insertcols!(gen_df, Symbol("block_min") => [minimum(r.m_surr) for r in res])
insertcols!(gen_df, Symbol("block_max") => [maximum(r.m_surr) for r in res])

# test curtailment - wavelet surrogates
res = @showprogress pmap(1:nrow(curtail_df)) do i
    surro_file = joinpath(wavelet_surros_dir, "curtail$(i)_surro_s$(scenario)_$(yr).txt")
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return curtail_independence(load_bus, i, scenario, yr, cfg, PrecomputedWLS(file=surro_file))
end

# save results
insertcols!(curtail_df, Symbol("wls_pvalue") => [r.pvalue for r in res])
insertcols!(curtail_df, Symbol("wls_est") => [r.m for r in res])
insertcols!(curtail_df, Symbol("wls_min") => [minimum(r.m_surr) for r in res])
insertcols!(curtail_df, Symbol("wls_max") => [maximum(r.m_surr) for r in res])

# test curtailment - block-shuffled surrogates
res = @showprogress pmap(1:nrow(curtail_df)) do i
    cfg = IndependenceTestConfig(k_embed=50, k_est=20, bin_width=0.2, nshuffles=200, w=23)

    return curtail_independence(load_bus, i, scenario, yr, cfg, BlockShuffle(2; shift=true))
end

# save results
insertcols!(curtail_df, Symbol("block_pvalue") => [r.pvalue for r in res])
insertcols!(curtail_df, Symbol("block_est") => [r.m for r in res])
insertcols!(curtail_df, Symbol("block_min") => [minimum(r.m_surr) for r in res])
insertcols!(curtail_df, Symbol("block_max") => [maximum(r.m_surr) for r in res])

# write results
CSV.write(joinpath(output_dir, "batt_res_s$(scenario)_$(yr)_$(parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])).csv"), batt_df)
CSV.write(joinpath(output_dir, "flow_res_s$(scenario)_$(yr)_$(parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])).csv"), flow_df)
CSV.write(joinpath(output_dir, "if_flow_res_s$(scenario)_$(yr)_$(parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])).csv"), if_flow_df)
CSV.write(joinpath(output_dir, "gen_res_s$(scenario)_$(yr)_$(parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])).csv"), gen_df)
CSV.write(joinpath(output_dir, "curtail_res_s$(scenario)_$(yr)_$(parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])).csv"), curtail_df)

rmprocs(workers())