# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

# load dependencies
using CairoMakie
using ColorSchemes: amp
using CSV
import Random

# include functions from other scripts
include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))
include(joinpath(@__DIR__, "..", "utils", "embedding.jl"))

# helper function
# fills param_df with optimal parameters, plots MSRE heatmaps and predicted time series
function compare_embedding_opts!(param_df::DataFrame, data_df::DataFrame, scaled_data_label::String, fig_opts)
    fig_ylabel, fig_title, fig_dir_prefix = fig_opts

    for id in 1:nrow(data_df)
        s = data_df[id, Symbol(scaled_data_label)]

        # try k = 25
        k = 25
        opt = EmbeddingOptions(max_d=8, max_τ=24, nbhd_type="mass", nbhd_param=k, w=23, metric=Cityblock())
        query_idx = collect(1+(opt.max_d-1)*opt.max_τ+opt.u : length(s)) # standardizes length of time series across all (d,τ)-combos

        mre25 = zeros(opt.max_τ, opt.max_d)
        for τ in 1:opt.max_τ, d in 1:opt.max_d
            for ct in 1:10
                local_pred = find_local_pred(s, query_idx, d, τ, opt)
                mre25[τ, d] += ragwitz_mre(s, query_idx, local_pred)
            end
        end
        mre25 ./= 10
        d_25 = findmin(mre25)[2][2]
        if d_25 == 1 # if embedding dimension 1, embedded "vector" is just value at previous timestep
            param_df[(param_df.id .== id) .& (param_df.k .== 25), :dim]   = [1]
            param_df[(param_df.id .== id) .& (param_df.k .== 25), :delay] = [-1]
        else
            param_df[(param_df.id .== id) .& (param_df.k .== 25), :dim]   = [d_25]
            param_df[(param_df.id .== id) .& (param_df.k .== 25), :delay] = [-findmin(mre25)[2][1]]
        end

        # try k = 50
        k = 50
        opt = EmbeddingOptions(max_d=8, max_τ=24, nbhd_type="mass", nbhd_param=k, w=23, metric=Cityblock())

        mre50 = zeros(opt.max_τ, opt.max_d)
        for τ in 1:opt.max_τ, d in 1:opt.max_d
            for ct in 1:10
                local_pred = find_local_pred(s, query_idx, d, τ, opt)
                mre50[τ, d] += ragwitz_mre(s, query_idx, local_pred)
            end
        end
        mre50 ./= 10
        d_50 = findmin(mre50)[2][2]
        if d_50 == 1 # if embedding dimension 1, embedded "vector" is just value at previous timestep
            param_df[(param_df.id .== id) .& (param_df.k .== 50), :dim]   = [1]
            param_df[(param_df.id .== id) .& (param_df.k .== 50), :delay] = [-1]
        else
            param_df[(param_df.id .== id) .& (param_df.k .== 50), :dim]   = [d_50]
            param_df[(param_df.id .== id) .& (param_df.k .== 50), :delay] = [-findmin(mre50)[2][1]]
        end
        
        min_err, max_err = extrema(vcat(mre25, mre50))
        if min_err != max_err
            fig = Figure(size=(1100, 450))
            g = fig[2,1] = GridLayout()
            ax21 = Axis(g[1,1], 
                        ylabel="dimension", 
                        ylabelsize=20, 
                        xlabel="delay", 
                        xlabelsize=20, 
                        yticks=1:opt.max_d, 
                        yticklabelsize=18, 
                        xticks=2:2:opt.max_τ, 
                        xticklabelsize=18, 
                        title="MRSE for local model predictions with k = 25", 
                        titlesize=20
            )
            colorrange = extrema(vcat(mre25, mre50))
            heatmap!(ax21, 1:opt.max_τ, 1:opt.max_d, mre25, colormap=:amp, colorrange=colorrange)
            scatter!(ax21, findmin(mre25)[2][1], findmin(mre25)[2][2], marker=:xcross, color=:black, markersize=:15)
            ax22 = Axis(g[1,2], 
                        ylabel="dimension", 
                        ylabelsize=20, 
                        xlabel="delay", 
                        xlabelsize=20, 
                        yticks=1:opt.max_d, 
                        yticklabelsize=18, 
                        xticks=2:2:opt.max_τ, 
                        xticklabelsize=18, 
                        title="MRSE for local model predictions with k = 50", 
                        titlesize=20
            )
            heatmap!(ax22, 1:opt.max_τ, 1:opt.max_d, mre50, colormap=:amp, colorrange=colorrange)
            scatter!(ax22, findmin(mre50)[2][1], findmin(mre50)[2][2], marker=:xcross, color=:black, markersize=:15)
            Colorbar(g[1,3], colormap=:amp, colorrange=colorrange, ticklabelsize=18)
            ax1 = Axis(fig[1,1],
                        ylabel=fig_ylabel,
                        xlabel="hour",
                        xgridvisible=false,
                        ygridvisible=false,
                        ylabelsize=20,
                        xlabelsize=20,
                        xticklabelsize=18,
                        yticklabelsize=18,
                        title=fig_title*" $(id)",
                        titlesize=20
            )
            lines!(ax1, s, linewidth=1.5, color=:gray22)
            rowsize!(fig.layout, 1, Auto(0.55))
            save(fig_dir_prefix*"$(id).png", fig)
        end
    end
end

# main script
# -----------

# base scenarios
scenario_yr_pairs = [(140, 2012), (290, 2002), (69, 2011)]
scenario, yr = scenario_yr_pairs[parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])]
Random.seed!(scenario+yr) # for reproducibility

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))

# output directories
output_dir = joinpath(@__DIR__, "embedding_params")
fig_dir = joinpath(@__DIR__, "embedding_params", "figures", "s$(scenario)_$(yr)")
if !isdir(fig_dir)
    mkdir(fig_dir)
end

# load shedding series
load_buses = [44]
ls_df = scale_load_shed_results(scenario, yr)[load_buses, :]
transform!(ls_df, :scaled_load_shed => ByRow(x -> x[JJA_idx]) => :scaled_load_shed)

n_load_buses = length(load_buses)
ls_param_df = DataFrame(scenario = scenario*ones(Int,2*n_load_buses), 
                        yr = yr*ones(Int,2*n_load_buses), 
                        id = repeat(1:n_load_buses, inner=2), 
                        k = repeat([25, 50], outer=n_load_buses), 
                        dim = zeros(Int,2*n_load_buses), 
                        delay = zeros(Int,2*n_load_buses)
)
fig_opts = (ylabel="load shed", title="Load bus", fig_dir_prefix=joinpath(fig_dir, "ls"))
compare_embedding_opts!(ls_param_df, ls_df, "scaled_load_shed", fig_opts)
ls_param_df.id = load_buses[ls_param_df.id]

# save output
if isfile(joinpath(output_dir, "ls_embed_params.csv"))
    CSV.write(joinpath(output_dir, "ls_embed_params.csv"), ls_param_df, append=true)
else
    CSV.write(joinpath(output_dir, "ls_embed_params.csv"), ls_param_df)
end

# battery state series
batt_df = scale_batt_results(scenario, yr)
transform!(batt_df, :scaled_state => ByRow(x -> x[JJA_idx]) => :scaled_state)
batt_param_df = DataFrame(scenario = scenario*ones(Int,2*nrow(batt_df)), 
                          yr = yr*ones(Int,2*nrow(batt_df)), 
                          id = repeat(1:nrow(batt_df), inner=2), 
                          k = repeat([25, 50], outer=nrow(batt_df)), 
                          dim = zeros(Int,2*nrow(batt_df)), 
                          delay = zeros(Int,2*nrow(batt_df))
)
fig_opts = (ylabel="state", title="Battery", fig_dir_prefix=joinpath(fig_dir, "batt"))
compare_embedding_opts!(batt_param_df, batt_df, "scaled_state", fig_opts)

# save output
if isfile(joinpath(output_dir, "batt_embed_params.csv"))
    CSV.write(joinpath(output_dir, "batt_embed_params.csv"), batt_param_df, append=true)
else
    CSV.write(joinpath(output_dir, "batt_embed_params.csv"), batt_param_df)
end

# flow series
branches = [99, 100] # last two branches are new HVDC lines
flow_df = scale_flow_results(scenario, yr)[branches, :]
transform!(flow_df, :scaled_flow => ByRow(x -> x[JJA_idx]) => :scaled_flow)
flow_param_df = DataFrame(scenario = scenario*ones(Int,2*nrow(flow_df)), 
                          yr = yr*ones(Int,2*nrow(flow_df)), 
                          id = repeat(1:nrow(flow_df), inner=2), 
                          k = repeat([25, 50], outer=nrow(flow_df)), 
                          dim = zeros(Int,2*nrow(flow_df)), 
                          delay = zeros(Int,2*nrow(flow_df))
)
fig_opts = (ylabel="flow", title="Branch", fig_dir_prefix=joinpath(fig_dir, "flow"))
compare_embedding_opts!(flow_param_df, flow_df, "scaled_flow", fig_opts)
flow_param_df.id = branches[flow_param_df.id]

# save output
if isfile(joinpath(output_dir, "flow_embed_params.csv"))
    CSV.write(joinpath(output_dir, "flow_embed_params.csv"), flow_param_df, append=true)
else
    CSV.write(joinpath(output_dir, "flow_embed_params.csv"), flow_param_df)
end

# IF flow series
if_flow_df = scale_if_flow_results(scenario, yr)
transform!(if_flow_df, :scaled_if_flow => ByRow(x -> x[JJA_idx]) => :scaled_if_flow)
if_flow_param_df = DataFrame(scenario = scenario*ones(Int,2*nrow(if_flow_df)), 
                             yr = yr*ones(Int,2*nrow(if_flow_df)), 
                             id = repeat(1:nrow(if_flow_df), inner=2), 
                             k = repeat([25, 50], outer=nrow(if_flow_df)), 
                             dim = zeros(Int,2*nrow(if_flow_df)), 
                             delay = zeros(Int,2*nrow(if_flow_df))
)
fig_opts = (ylabel="IF flow", title="Interface", fig_dir_prefix=joinpath(fig_dir, "if"))
compare_embedding_opts!(if_flow_param_df, if_flow_df, "scaled_if_flow", fig_opts)

# save output
if isfile(joinpath(output_dir, "if_flow_embed_params.csv"))
    CSV.write(joinpath(output_dir, "if_flow_embed_params.csv"), if_flow_param_df, append=true)
else
    CSV.write(joinpath(output_dir, "if_flow_embed_params.csv"), if_flow_param_df)
end

# gen series
gens = 1:29 # dispatchable energy sources (hydro & imported)
gen_df = scale_gen_results(scenario, yr)[gens, :]
transform!(gen_df, :scaled_gen => ByRow(x -> x[JJA_idx]) => :scaled_gen)
gen_param_df = DataFrame(scenario = scenario*ones(Int,2*nrow(gen_df)), 
                         yr = yr*ones(Int,2*nrow(gen_df)), 
                         id = repeat(1:nrow(gen_df), inner=2), 
                         k = repeat([25, 50], outer=nrow(gen_df)), 
                         dim = zeros(Int,2*nrow(gen_df)), 
                         delay = zeros(Int,2*nrow(gen_df))
)
fig_opts = (ylabel="generation", title="Generator", fig_dir_prefix=joinpath(fig_dir, "gen"))
compare_embedding_opts!(gen_param_df, gen_df, "scaled_gen", fig_opts)

# save output
if isfile(joinpath(output_dir, "gen_embed_params.csv"))
    CSV.write(joinpath(output_dir, "gen_embed_params.csv"), gen_param_df, append=true)
else
    CSV.write(joinpath(output_dir, "gen_embed_params.csv"), gen_param_df)
end

# curtailment series - non-dispatchable energy sources
curtail_df = scale_curtail_results(scenario, yr)
transform!(curtail_df, :scaled_curtail => ByRow(x -> x[JJA_idx]) => :scaled_curtail)
curtail_param_df = DataFrame(scenario = scenario*ones(Int,2*nrow(curtail_df)), 
                             yr = yr*ones(Int,2*nrow(curtail_df)), 
                             id = repeat(1:nrow(curtail_df), inner=2), 
                             k = repeat([25, 50], outer=nrow(curtail_df)), 
                             dim = zeros(Int,2*nrow(curtail_df)), 
                             delay = zeros(Int,2*nrow(curtail_df))
)
fig_opts = (ylabel="curtailment", title="Curtailed generator", fig_dir_prefix=joinpath(fig_dir, "curtail"))
compare_embedding_opts!(curtail_param_df, curtail_df, "scaled_curtail", fig_opts)

# save output
if isfile(joinpath(output_dir, "curtail_embed_params.csv"))
    CSV.write(joinpath(output_dir, "curtail_embed_params.csv"), curtail_param_df, append=true)
else
    CSV.write(joinpath(output_dir, "curtail_embed_params.csv"), curtail_param_df)
end