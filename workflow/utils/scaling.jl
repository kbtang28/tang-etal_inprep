using DataFrames
using DelimitedFiles

import CSV

include(joinpath(@__DIR__, "other.jl"))

acorn_dir = joinpath(@__DIR__, "..", "..", "ACORN")
acorn_results_dir = joinpath(acorn_dir, "results")

# scales relevant generation results (hydro, imported, wind, and solar)
function scale_gen_results(scenario::Int, yr::Int)
    # get relevant DU factors
    _, _, _, wind_cap_sf, solar_cap_sf, _, du_scenario = get_du_factors(scenario)
    
    # read in relevant generator info
    gen_df = DataFrame(CSV.File(joinpath(acorn_dir, "gen_info.csv")))
    transform!(gen_df, 
               :max_gen => ByRow(max_gen -> max_gen == 0.0 ? Vector{Float64}(undef, 8760) : max_gen) => :max_gen
    )

    # scale MS hydro max gen capacity
    mshydro_df = DataFrame(CSV.File(joinpath(acorn_dir, "hydro", "nypaMosesSaundersEnergy.climate.change.csv"); 
                                    header=true, skipto=(yr-1953)*48+2, footerskip=(2019-yr)*48)
    )
    ms_cap = mshydro_df[:, Symbol("nypaMosesSaundersEnergy.$(du_scenario)")]

    qm_to_numdays = readdlm(joinpath(@__DIR__, "qm_to_num_days.csv"), ',', Int; skipstart=1)
    qm_to_nhours = qm_to_numdays[:, 2]*24

    ms_cap_sf = maximum(ms_cap ./ qm_to_nhours / gen_df.max_gen[5])
    if ms_cap_sf > 1 # increase max gen capacity of Moses Saunders gen if necessary
        gen_df[gen_df.id .== 5, :max_gen] *= ms_cap_sf
    end

    # read in gen results
    gen_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "gen_$(yr).csv"), ',', Float64)
    gen_df.gen = [gen_array[i, :] for i in 1:size(gen_array, 1)]

    # get wind max gen
    wind = readdlm(joinpath(acorn_dir, "wind", "wind_$(yr).csv"), ',', Float64)[:, 2:end]
    wind_cap = wind_cap_sf * round.(wind; digits=2)
    wind_cap[wind_cap .== 0.0] .= 1.0 # to avoid division by 0

    gen_df[gen_df.type .== "Wind", :max_gen] = [wind_cap[i, :] for i in 1:size(wind_cap, 1)]

    # get solar max gen
    solar = readdlm(joinpath(acorn_dir, "solar", "Scenario$(du_scenario)", "solarUPV$(yr).csv"), ',', Float64)[:, 2:end]
    solar_cap = solar_cap_sf * round.(solar; digits=2)
    solar_cap[solar_cap .== 0.0] .= 1.0

    gen_df[gen_df.type .== "SolarUPV", :max_gen] = [solar_cap[i, :] for i in 1:size(solar_cap, 1)]

    # scale gen
    select!(gen_df, :id, :bus_intern, :bus_extern, :type, :gen, [:gen, :max_gen] => ByRow((gen, max_gen) -> gen ./ max_gen) => :scaled_gen)
    
    # keep only relevant gen results and cols
    subset!(gen_df, :type => type -> type .∈ (["Hydro", "Import", "Wind", "SolarUPV"], ))

    return gen_df
end

# scales wind and solar curtailment results
function scale_curtail_results(scenario::Int, yr::Int)
    # get relevant DU factors
    _, _, _, wind_cap_sf, solar_cap_sf, _, du_scenario = get_du_factors(scenario)

    # read in relevant generator info
    gen_df = DataFrame(CSV.File(joinpath(acorn_dir, "gen_info.csv")))
    curtail_df = gen_df[(gen_df.type .== "Wind") .|| (gen_df.type .== "SolarUPV"), :]
    transform!(curtail_df,
               :max_gen => ByRow(max_gen -> max_gen == 0.0 ? Vector{Float64}(undef, 8760) : max_gen) => :max_gen
    )

    # read in curtailment results
    wc_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "wc_$(yr).csv"), ',', Float64)
    sc_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "sc_$(yr).csv"), ',', Float64)
    curtail_df.curtail = vcat([wc_array[i, :] for i in 1:size(wc_array, 1)], [sc_array[i, :] for i in 1:size(sc_array, 1)])

    # get wind max gen
    wind = readdlm(joinpath(acorn_dir, "wind", "wind_$(yr).csv"), ',', Float64)[:, 2:end]
    wind_cap = wind_cap_sf * round.(wind; digits=2)
    wind_cap[wind_cap .== 0.0] .= 1.0 # to avoid division by 0

    curtail_df[curtail_df.type .== "Wind", :max_gen] = [wind_cap[i, :] for i in 1:size(wind_cap, 1)]

    # get solar max gen
    solar = readdlm(joinpath(acorn_dir, "solar", "Scenario$(du_scenario)", "solarUPV$(yr).csv"), ',', Float64)[:, 2:end]
    solar_cap = solar_cap_sf * round.(solar; digits=2)
    solar_cap[solar_cap .== 0.0] .= 1.0

    curtail_df[curtail_df.type .== "SolarUPV", :max_gen] = [solar_cap[i, :] for i in 1:size(solar_cap, 1)]

    # scale curtailment
    select!(curtail_df, :id, :bus_intern, :bus_extern, :type, :curtail, [:curtail, :max_gen] => ByRow((curtail, max_gen) -> curtail ./ max_gen) =>:scaled_curtail)

    return curtail_df
end

# scales battery storage results
function scale_batt_results(scenario::Int, yr::Int)
    # get relevant DU factors
    _, _, _, _, _, batt_cap_sf, _ = get_du_factors(scenario)

    # read in battery info
    batt_df = DataFrame(CSV.File(joinpath(acorn_dir, "storage_info.csv")))
    batt_df.storage_cap = (batt_df.capacity .* batt_df.duration) * batt_cap_sf

    # scale battery state results
    state_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "battstate_$(yr).csv"), ',', Float64)[:, 2:end]
    batt_df.state = [state_array[i, :] for i in 1:size(state_array, 1)]
    batt_df.scaled_state = batt_df.state ./ batt_df.storage_cap
    
    return batt_df
end

# scales branch flow results
function scale_flow_results(scenario, yr)
    # get relevant DU factors
    _, _, _, _, _, _, du_scenario = get_du_factors(scenario)

    # read in regular line info
    flow_df = DataFrame(CSV.File(joinpath(acorn_dir, "branches.csv")))
    flow_df.dc = fill(false, nrow(flow_df))

    # read in regular line flow results
    flow_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "flow_$(yr).csv"), ',', Float64)
    flow_df.flow = [flow_array[i, :] for i in 1:size(flow_array, 1)]

    # scale regular line flow results - unconstrained branches
    scaled_flow_array = zeros(Float64, size(flow_array))    
    unconst_branch_ids = vcat(1:3, 5:7, 9:57, 59, 60, 62:67, 69, 71:94)
    scaled_flow_array[unconst_branch_ids, :] = flow_array[unconst_branch_ids, :] ./ maximum(abs.(flow_array[unconst_branch_ids, :]), dims=2)

    # scale regular line flow results - constrained branches
    rated_branch_ids = [58, 61, 70]
    scaled_flow_array[rated_branch_ids, :] = flow_array[rated_branch_ids, :] ./ flow_df.rating[rated_branch_ids]

    if_branch_ids = [8; 4; 68]
    if_ids = [6; 7; 9]
    iflim_dn_array = readdlm(joinpath(acorn_dir, "if_limits", "iflimdn_$(yr)_$(du_scenario).csv"), ',', Float64)
    iflim_up_array = readdlm(joinpath(acorn_dir, "if_limits", "iflimup_$(yr)_$(du_scenario).csv"), ',', Float64)
    for (if_branch_id, if_id) in zip(if_branch_ids, if_ids)
        pos_flow_t = view(flow_array, if_branch_id, :) .>= 0.0
        scaled_flow_array[if_branch_id, pos_flow_t] = flow_array[if_branch_id, pos_flow_t] ./ iflim_up_array[if_id, pos_flow_t]
        scaled_flow_array[if_branch_id, .!pos_flow_t] = -flow_array[if_branch_id, .!pos_flow_t] ./ iflim_dn_array[if_id, .!pos_flow_t]
    end

    flow_df.scaled_flow = [scaled_flow_array[i, :] for i in 1:size(scaled_flow_array, 1)]

    # read in DC line info
    dc_line_df = DataFrame(CSV.File(joinpath(acorn_dir, "dclines.csv")))
    dc_line_df.id = flow_df.id[end]+1:flow_df.id[end]+nrow(dc_line_df)
    dc_line_df.dc = fill(true, nrow(dc_line_df))
    
    # read in DC line results
    gen_df = DataFrame(CSV.File(joinpath(acorn_dir, "gen_info.csv"); skipto=69))
    gen_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "gen_$(yr).csv"), ',', Float64; skipstart=67)
    gen_array = gen_array[gen_df.type .== "tDummy", :]

    # scale DC line results
    dc_line_df.flow = [gen_array[i, :] for i in 1:size(gen_array, 1)]
    transform!(dc_line_df,
               [:rating, :flow] => ByRow((rating, flow) -> flow/rating) => :scaled_flow
    )

    # concatenate
    flow_df = vcat(flow_df, dc_line_df)

    return flow_df
end

# scales IF flow results
function scale_if_flow_results(scenario, yr)
    # get relevant DU factors
    _, _, _, _, _, _, du_scenario = get_du_factors(scenario)

    if_flow_df = DataFrame(:id => 1:15)

    # read in IF flow results
    if_flow_array = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "ifsum_$(yr).csv"), ',', Float64)
    if_flow_df.if_flow = [if_flow_array[i, :] for i in 1:size(if_flow_array, 1)]

    # read in IF limits
    iflim_dn_array = readdlm(joinpath(acorn_dir, "if_limits", "iflimdn_$(yr)_$(du_scenario).csv"), ',', Float64)
    iflim_up_array = readdlm(joinpath(acorn_dir, "if_limits", "iflimup_$(yr)_$(du_scenario).csv"), ',', Float64)

    scaled_if_flow = Vector{Float64}[]
    for i in 1:15
        pos_flow_t = (view(if_flow_array, i, :) .>= 0.0)
        if_flow_array[i, pos_flow_t] = if_flow_array[i, pos_flow_t] ./ iflim_up_array[i, pos_flow_t]
        if_flow_array[i, .!pos_flow_t] = -if_flow_array[i, .!pos_flow_t] ./ iflim_dn_array[i, .!pos_flow_t]

        push!(scaled_if_flow, if_flow_array[i, :])
    end
    if_flow_df.scaled_if_flow = scaled_if_flow

    return if_flow_df
end

# calculates load at each bus using node balance constraint
function calc_load(scenario, yr)
    # read in load shedding results
    load_shed = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "loadshed_$(yr).csv"), ',', Float64)
    
    # read in branch info and flow results
    flow_df = DataFrame(CSV.File(joinpath(acorn_dir, "branches.csv")))
    flow = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "flow_$(yr).csv"), ',', Float64)

    # read in generator info and generation results
    gen_df = DataFrame(CSV.File(joinpath(acorn_dir, "gen_info.csv")))
    gen = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "gen_$(yr).csv"), ',', Float64)
    
    # read in storage info and results
    storage_df = DataFrame(CSV.File(joinpath(acorn_dir, "storage_info.csv")))
    charge = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "charege_$(yr).csv"), ',', Float64)
    discharge = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "disch_$(yr).csv"), ',', Float64)

    # process info and calculate load for each bus
    load = zeros(size(load_shed))
    
    n_bus = size(load_shed, 1)
    for b = 1:n_bus
        fbranch_id = (flow_df.fbus_intern .== b)
        tbranch_id = (flow_df.tbus_intern .== b)
        gen_id = (gen_df.bus_intern .== b)
        batt_id = (storage_df.bus_intern .== b)

        load[b, :] = (- vec(sum(view(flow, fbranch_id, :), dims=1))     # power flowing out via branches 
                      + vec(sum(view(flow, tbranch_id, :), dims=1))     # power flowing in via branches
                      + vec(sum(view(gen, gen_id, :), dims=1))          # power generated at bus
                      + vec(sum(view(discharge, batt_id, :), dims=1))   # power discharged from battery at bus
                      - vec(sum(view(charge, batt_id, :), dims=1))      # power used to charge battery at bus
                      + view(load_shed, b, :)                           # load shed at bus
        )
    end
    load[abs.(load) .< 1e-4] .= 0.0

    return load
end

# ... same as above but avoids reading load shedding results twice
function _calc_load(scenario, yr, load_shed)
    # read in branch info and flow results
    flow_df = DataFrame(CSV.File(joinpath(acorn_dir, "branches.csv")))
    flow = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "flow_$(yr).csv"), ',', Float64)

    # read in generator info and generation results
    gen_df = DataFrame(CSV.File(joinpath(acorn_dir, "gen_info.csv")))
    gen = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "gen_$(yr).csv"), ',', Float64)
    
    # read in storage info and results
    storage_df = DataFrame(CSV.File(joinpath(acorn_dir, "storage_info.csv")))
    charge = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "charege_$(yr).csv"), ',', Float64)
    discharge = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "disch_$(yr).csv"), ',', Float64)

    # process info and calculate load for each bus
    load = zeros(size(load_shed))
    
    n_bus = size(load_shed, 1)
    for b = 1:n_bus
        fbranch_id = (flow_df.fbus_intern .== b)
        tbranch_id = (flow_df.tbus_intern .== b)
        gen_id = (gen_df.bus_intern .== b)
        batt_id = (storage_df.bus_intern .== b)

        load[b, :] = (- vec(sum(view(flow, fbranch_id, :), dims=1))     # power flowing out via branches 
                      + vec(sum(view(flow, tbranch_id, :), dims=1))     # power flowing in via branches
                      + vec(sum(view(gen, gen_id, :), dims=1))          # power generated at bus
                      + vec(sum(view(discharge, batt_id, :), dims=1))   # power discharged from battery at bus
                      - vec(sum(view(charge, batt_id, :), dims=1))      # power used to charge battery at bus
                      + view(load_shed, b, :)                           # load shed at bus
        )
    end
    load[abs.(load) .< 1e-4] .= 0.0

    return load
end

# scales load shedding
function scale_load_shed_results(scenario, yr)
    # read in load shedding results
    load_shed = readdlm(joinpath(acorn_results_dir, "Scenario$(scenario)", "loadshed_$(yr).csv"), ',', Float64)
    
    # get total load
    load = _calc_load(scenario, yr, load_shed)

    # scale load shedding
    scaled_load_shed = zeros(size(load_shed))
    scaled_load_shed[load .> 0.0] = load_shed[load .> 0.0] ./ load[load .> 0.0]

    # read in bus info and add scaled load shedding
    load_shed_df = DataFrame(CSV.File(joinpath(acorn_dir, "bus_ids.csv")))
    load_shed_df.load_shed = [load_shed[i, :] for i in 1:size(load_shed, 1)]
    load_shed_df.scaled_load_shed = [scaled_load_shed[i, :] for i in 1:size(scaled_load_shed, 1)]

    return load_shed_df
end