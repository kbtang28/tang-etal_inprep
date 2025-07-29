using CSV
using DataFrames
using DelimitedFiles

# returns DU parameter combinations for all scenarios
function get_du_factors_df()::DataFrame
    header = [:temp_increase; 
              :building_elec_rate; 
              :ev_elec_rate; 
              :wind_cap_sf; 
              :solar_cap_sf; 
              :batt_cap_sf; 
              :du_scenario]
    df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "ACORN", "du_factors.csv"); header=false), header)
    df[!, :du_scenario] = convert.(Int, df[!, :du_scenario])
    sort!(df, :du_scenario)

    return df
end

# returns DU parameter combination for specific scenario
function get_du_factors(scenario::Int64)
    df = get_du_factors_df()

    return df[scenario, :temp_increase], df[scenario, :building_elec_rate], df[scenario, :ev_elec_rate], df[scenario, :wind_cap_sf], df[scenario, :solar_cap_sf], df[scenario, :batt_cap_sf], df[scenario, :du_scenario]
end

# converts internal bus ID to external bus ID
function get_extern_id(intern_id)
    bus_ids = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "ACORN", "bus_ids.csv")))

    return bus_ids[bus_ids.intern_id .== intern_id, :extern_id][1]
end

# converts external bus ID to internal bus ID
function get_intern_id(extern_id)
    bus_ids = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "ACORN", "bus_ids.csv")))

    return bus_ids[bus_ids.extern_id .== extern_id, :intern_id][1]
end

# returns the time indices for given month
function get_month_idx(month)
    qm_to_num_days = readdlm(joinpath(@__DIR__, "qm_to_num_days.csv"), ',', Int; skipstart=1)
    
    if month == 1
        days = sum(qm_to_num_days[1:4, 2])
        return collect(1:days*24)
    else
        qm_start = (month - 1)*4 + 1
        days = sum(qm_to_num_days[qm_start:qm_start+3, 2])
        t0 = sum(qm_to_num_days[1:qm_start-1, 2])*24
        return collect(t0+1:t0+days*24)
    end    
end

function get_qm_idx(qm)
    qm_to_num_days = readdlm(joinpath(@__DIR__, "qm_to_num_days.csv"), ',', Int; skipstart=1)

    if qm == 1
        days = qm_to_num_days[1, 2]
        return collect(1:days*24)
    else
        days = qm_to_num_days[qm, 2]
        t0 = sum(qm_to_num_days[1:(qm-1), 2])*24
        return collect(t0+1:t0+days*24)
    end
end