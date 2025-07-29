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
    using DelimitedFiles
    using Random
    using TimeseriesSurrogates

    include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))

    Random.seed!(myid()) # for reproducibility
end

# main script
# -----------

scenario_yr_pairs = [(140, 2012), (290, 2002), (69, 2011)]
scenario, yr = scenario_yr_pairs[parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])]
output_dir = joinpath(@__DIR__, "precomputed_wls_surrogates")

nsurro = 200 # number of surrogates to generate
method = WLS(IAAFT(), rescale=true)

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))

batt_df = scale_batt_results(scenario, yr)
transform!(batt_df, :scaled_state => ByRow(x -> x[JJA_idx]) => :scaled_state)
for i in 1:nrow(batt_df)
    s = batt_df[i, :scaled_state]
    ŝ = pmap(1:nsurro) do j
        surrogate(s, method)
    end
        
    writedlm(joinpath(output_dir, "batt$(i)_surro_s$(scenario)_$(yr).txt"), stack(ŝ))
end

branches = [99, 100] # last two branches are new HVDC lines
flow_df = scale_flow_results(scenario, yr)[branches, :]
transform!(flow_df, :scaled_flow => ByRow(x -> x[JJA_idx]) => :scaled_flow)
for (i, id) in enumerate(flow_df.id)
    s = flow_df[i, :scaled_flow]
    s = s ./ 2 # scale to have theoretical range ≤ 1

    ŝ = pmap(1:nsurro) do j
        surrogate(s, method)
    end

    writedlm(joinpath(output_dir, "flow$(id)_surro_s$(scenario)_$(yr).txt"), stack(ŝ))
end

if_flow_df = scale_if_flow_results(scenario, yr)
transform!(if_flow_df, :scaled_if_flow => ByRow(x -> x[JJA_idx]) => :scaled_if_flow)
for i in 1:nrow(if_flow_df)
    s = if_flow_df[i, :scaled_if_flow]
    s = s ./ 2 # scale to have theoretical range ≤ 1

    ŝ = pmap(1:nsurro) do j
        surrogate(s, method)
    end

    writedlm(joinpath(output_dir, "if_flow$(i)_surro_s$(scenario)_$(yr).txt"), stack(ŝ))
end

gens = 1:29 # dispatchable energy sources (hydro & imported)
gen_df = scale_gen_results(scenario, yr)[gens, :]
transform!(gen_df, :scaled_gen => ByRow(x -> x[JJA_idx]) => :scaled_gen)
for i in 1:nrow(gen_df)
    s = gen_df[i, :scaled_gen]
    ŝ = pmap(1:nsurro) do j
        surrogate(s, method)
    end

    writedlm(joinpath(output_dir, "gen$(i)_surro_s$(scenario)_$(yr).txt"), stack(ŝ))
end

curtail_df = scale_curtail_results(scenario, yr)
transform!(curtail_df, :scaled_curtail => ByRow(x -> x[JJA_idx]) => :scaled_curtail)
for i in 1:nrow(curtail_df)
    s = curtail_df[i, :scaled_curtail]
    ŝ = pmap(1:nsurro) do j
        surrogate(s, method)
    end

    writedlm(joinpath(output_dir, "curtail$(i)_surro_s$(scenario)_$(yr).txt"), stack(ŝ))
end

rmprocs(workers())