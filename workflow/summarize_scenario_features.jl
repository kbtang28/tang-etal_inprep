# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# load dependencies and helper functions
using Associations
using DataFrames
using DelimitedFiles
using ProgressMeter

import StatsBase: mean

include(joinpath(@__DIR__, "utils", "scaling.jl"))
include(joinpath(@__DIR__, "utils", "embedding.jl"))

# utility function to estimate conditional entropy H(Y⁺ | Y) - theoretical upper bound on TE(X → Y)
function estimate_ls_conditional_entropy(s)
    # optimize embedding for series
    opt = EmbeddingOptions(max_d=8, max_τ=24, nbhd_type="mass", nbhd_param=50, w=23)
    d, τ = optimize_embedding(s, opt; n_itr=20)

    # truncate for binning
    s[s .< 0.0] .= 0.0
    s[s .> 0.99] .= 0.99

    # embed series
    X = genembed(s, collect(0:τ:(d-1)*τ))[1:end-1]
    Xf = s[( 2+(d-1)*(-τ) ):end]

    # discretize and estimate conditional entropy of Xf given X
    precise = true
    encodingXf = RectangularBinEncoding(FixedRectangularBinning(0.0:0.2:1.0, 1, precise))
    encodingX = RectangularBinEncoding(FixedRectangularBinning(0.0:0.2:1.0, dimension(X), precise))
    disc = CodifyPoints(encodingXf, encodingX)
    est = JointProbabilities(ConditionalEntropyShannon(), disc)
    return association(est, Xf, X)
end

JJA_idx = vcat(get_month_idx(6), get_month_idx(7), get_month_idx(8))

acorn_dir = joinpath(@__DIR__, "..", "ACORN")
output_dir = joinpath(@__DIR__, "..", "output")
VML_data_dir = joinpath(@__DIR__, "..", "..", "..", "..", "shared", "vs498_0001", "Vivienne_SharedData") # full dataset from Vivienne

# main script
# -----------

df = DataFrame(:scenario => repeat(1:300, outer=22),
               :yr => repeat(1998:2019, inner=300),
               :avg_temp_anomaly => Vector{Float64}(undef, 6600), # relative to historical baseline
               :avg_wind_anomaly => Vector{Float64}(undef, 6600), # relative to historical baseline
               :avg_solar_anomaly => Vector{Float64}(undef, 6600), # relative to historical baseline
               :avg_hydro_anomaly => Vector{Float64}(undef, 6600), # relative to historical baseline
               :avg_load_anomaly => Vector{Float64}(undef, 6600), # relative to historical baseline
               :prop_demand_unmet => Vector{Float64}(undef, 6600), # at target bus
               :prop_demand_unmet_hrs => Vector{Float64}(undef, 6600), # at target bus
               :prop_curtailed => Vector{Float64}(undef, 6600), # prop. wind and solar energy curtailed
               :prop_congested_pos_hrs => Vector{Float64}(undef, 6600), # E → G congestion
               :prop_congested_neg_hrs => Vector{Float64}(undef, 6600), # G → E congestion
               :prop_congested_pos_with_unmet_demand_hrs => Vector{Float64}(undef, 6600),
               :conditional_entropy => Vector{Float64}(undef, 6600), # of power shortage indicator at target bus
               :prop_zoneJK_with_zoneG_shortage => Vector{Float64}(undef, 6600)
)

# add battery capacity scaling factor
du_factors_df = get_du_factors_df()
insertcols!(du_factors_df, :scenario => 1:300)
leftjoin!(df, du_factors_df[:, [:scenario, :batt_cap_sf]], on=:scenario)
dropmissing!(df) # none missing but gets rid of type error

# baseline temp
temp_dir = joinpath(acorn_dir, "temp")
temp_dict = Dict{Int, Vector{Float64}}()
for yr in 1998:2019
    temp_series = zeros(Float64, length(JJA_idx))
    for zone in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
        temp_series .+= vec(readdlm(joinpath(temp_dir, "temp_zone$(zone)_$(yr).txt"), skipstart=1))[JJA_idx]
    end
    temp_dict[yr] = temp_series/11 # averaged over zones
end
baseline_temp = sum(values(temp_dict))./22 # averaged over 22-yr sim period
writedlm(joinpath(output_dir, "baseline", "baseline_temp.txt"), baseline_temp)

temp_anomaly_dict = Dict{Int, Float64}() # store since temp increase modeled as step change
for yr in 1998:2019
    temp_anomaly_dict[yr] = mean(temp_dict[yr] .- baseline_temp)
end

# baseline wind availability
wind_dir = joinpath(acorn_dir, "wind")
wind_dict = Dict{Int, Vector{Float64}}()
for yr in 1998:2019
    wind_series = readdlm(joinpath(wind_dir, "wind_$(yr).csv"), ',')
    wind_series = round.(wind_series[:, JJA_idx.+1], digits=2) # first col is bus ID
    wind_dict[yr] = vec(sum(wind_series, dims=1)) # total wind capacity across all buses
end
baseline_wind = sum(values(wind_dict))./22 # averaged over 22-yr sim period
writedlm(joinpath(output_dir, "baseline", "baseline_wind.txt"), baseline_wind)

# baseline solar availability
solar_dir = joinpath(acorn_dir, "solar")
baseline_solar = zeros(length(JJA_idx))
for yr in 1998:2019
    solar_series_upv = readdlm(joinpath(solar_dir, "Scenario0", "solarUPV$(yr)_0.csv"), ',')
    solar_series_upv = round.(solar_series_upv[:, JJA_idx.+1], digits=2) # first col is bus ID
    baseline_solar .+= vec(sum(solar_series_upv, dims=1))
end
baseline_solar = baseline_solar./22 # averaged over 22-yr sim period
writedlm(joinpath(output_dir, "baseline", "baseline_solar.txt"), baseline_solar)

solar_dict = Dict{Int, Dict{Int, Vector{Float64}}}() # first key is DU scenario, second key is yr
for du_scenario in unique(du_factors_df.du_scenario)
    dict_tmp = Dict{Int, Vector{Float64}}()
    for yr in 1998:2019
        solar_series_upv = readdlm(joinpath(solar_dir, "Scenario$(du_scenario)", "solarUPV$(yr).csv"), ',')
        solar_series_upv = round.(solar_series_upv[:, JJA_idx.+1], digits=2) # first col is bus ID
        dict_tmp[yr] = vec(sum(solar_series_upv, dims=1))
    end
    solar_dict[du_scenario] = dict_tmp
end

# baseline hydro availability
hydro_dir = joinpath(acorn_dir, "hydro")
MS_hydro_df = DataFrame(CSV.File(joinpath(hydro_dir, "nypaMosesSaundersEnergy.climate.change.csv")))
N_hydro_df = DataFrame(CSV.File(joinpath(hydro_dir, "nypaNiagaraEnergy.climate.change.csv")))
baseline_hydro = zeros(12) # quarter-monthly series
for yr in 1998:2019
    # sum across two hydro generators
    baseline_hydro .+= MS_hydro_df[(MS_hydro_df.Year .== yr) .&& (in([6, 7, 8]).(MS_hydro_df.Month)), :nypaMosesSaundersEnergy]
    baseline_hydro .+= N_hydro_df[(N_hydro_df.Year .== yr) .&& (in([6, 7, 8]).(N_hydro_df.Month)), :nypaNiagaraEnergy]
end
baseline_hydro = baseline_hydro./22 # averaged over 22-yr sim period
writedlm(joinpath(output_dir, "baseline", "baseline_hydro.txt"), baseline_hydro)

hydro_anomaly_dict = Dict{Int, Dict{Int, Float64}}() # first key is DU scenario, second key is yr
for du_scenario in unique(du_factors_df.du_scenario)
    dict_tmp = Dict{Int, Float64}()
    for yr in 1998:2019
        hydro_series_MS = MS_hydro_df[(MS_hydro_df.Year .== yr) .&& (in([6, 7, 8]).(MS_hydro_df.Month)), Symbol("nypaMosesSaundersEnergy.$(du_scenario)")]
        hydro_series_N = N_hydro_df[(N_hydro_df.Year .== yr) .&& (in([6, 7, 8]).(N_hydro_df.Month)), Symbol("nypaNiagaraEnergy.$(du_scenario)")]
        dict_tmp[yr] = mean((hydro_series_MS .+ hydro_series_N) .- baseline_hydro)
    end
    hydro_anomaly_dict[du_scenario] = dict_tmp
end

# baseline load
bldg_elec_rate = 0.92 # from Vivienne's mainbasecaller.m
ev_elec_rate = 0.9 # from Vivienne's mainbasecaller.m

baseline_load = zeros(length(JJA_idx))
for yr in 1998:2019
    base_load = readdlm(joinpath(VML_data_dir, "BaseLoad", "Scenario0", "simload_$(yr).csv"), ',')
    base_load = vec(sum(base_load[:, JJA_idx], dims=1))

    res_load = readdlm(joinpath(VML_data_dir, "ResLoad", "Scenario0", "ResLoad_Bus_$(yr).csv"), ',')
    res_load = vec(sum(res_load[:, JJA_idx.+1], dims=1)) # first col is bus ID

    com_load = readdlm(joinpath(VML_data_dir, "ComLoad", "Scenario0", "ComLoad_Bus_$(yr).csv"), ',')
    com_load = vec(sum(com_load[:, JJA_idx.+1], dims=1)) # first col is bus ID

    ev_load = readdlm(joinpath(acorn_dir, "ev_load.csv"), ',')
    ev_load = vec(sum(ev_load[:, JJA_idx.+1], dims=1)) # first col is bus ID

    btm_solar = readdlm(joinpath(solar_dir, "Scenario0", "solarDPV$(yr)_0.csv"), ',')
    btm_solar = vec(sum(btm_solar[:, JJA_idx.+1], dims=1)) # first col is bus ID

    small_hydro = readdlm(joinpath(hydro_dir, "smallhydrogen.csv"), ',')
    small_hydro = vec(sum(small_hydro[:, JJA_idx], dims=1))

    baseline_load .+= base_load .+ bldg_elec_rate.*(res_load .+ com_load) .+ ev_elec_rate.*ev_load .- btm_solar .- small_hydro
end
baseline_load = baseline_load./22 # averaged over 22-yr sim period
writedlm(joinpath(output_dir, "baseline", "baseline_load.txt"), baseline_load)

@showprogress for scenario in 1:300
    temp_increase, _, _, wind_cap_sf, solar_cap_sf, _, du_scenario = get_du_factors(scenario)
    for yr in 1998:2019
        idx = findall((df.scenario .== scenario) .&& (df.yr .== yr))[1] # will only be one match

        # anomaly in temp, renewable availability, and load relative to baseline
        df[idx, :avg_temp_anomaly] = temp_anomaly_dict[yr] + temp_increase
        df[idx, :avg_wind_anomaly] = mean(wind_cap_sf.*wind_dict[yr] .- baseline_wind)
        df[idx, :avg_solar_anomaly] = mean(solar_cap_sf.*solar_dict[du_scenario][yr] .- baseline_solar)
        df[idx, :avg_hydro_anomaly] = hydro_anomaly_dict[du_scenario][yr]
        
        load_mat = calc_load(scenario, yr)[:, JJA_idx]
        df[idx, :avg_load_anomaly] = mean(vec(sum(load_mat, dims=1)) .- baseline_load)

        # unmet demand at target bus
        ls_mat = readdlm(joinpath(acorn_dir, "results", "Scenario$(scenario)", "loadshed_$(yr).csv"), ',', Float64)[:, JJA_idx]
        ls_hrs = (view(ls_mat, 44, :) .> 0.001) # hrs with unmet demand
        df[idx, :prop_demand_unmet] = round(sum(view(ls_mat, 44, :)) / sum(view(load_mat, 44, :)); digits=4)
        df[idx, :prop_demand_unmet_hrs] = sum(ls_hrs)/length(JJA_idx)

        # wind and solar curtailment
        sc_mat = readdlm(joinpath(acorn_dir, "results", "Scenario$(scenario)", "sc_$(yr).csv"), ',', Float64)[:, JJA_idx]
        wc_mat = readdlm(joinpath(acorn_dir, "results", "Scenario$(scenario)", "wc_$(yr).csv"), ',', Float64)[:, JJA_idx]
        wind_and_solar_total = sum(solar_cap_sf.*solar_dict[du_scenario][yr]) + sum(wind_cap_sf.*wind_dict[yr])
        df[idx, :prop_curtailed] = round(( sum(sc_mat) + sum(wc_mat) )/wind_and_solar_total; digits=4)

        # E-G IF congestion
        if6_flow = scale_if_flow_results(scenario, yr)[6, :scaled_if_flow][JJA_idx]
        congested_pos_hrs = (if6_flow .>= 0.999) # hrs with congestion in E → G direction
        df[idx, :prop_congested_pos_hrs] = sum(congested_pos_hrs)/length(JJA_idx)
        df[idx, :prop_congested_neg_hrs] = sum(if6_flow .<= -0.999)/length(JJA_idx)
        df[idx, :prop_congested_pos_with_unmet_demand_hrs] = sum(congested_pos_hrs .&& ls_hrs)/sum(congested_pos_hrs)
        
        # co-occurence of load shedding in zones J/K and zone G
        ls_df = scale_load_shed_results(scenario, yr)
        ls44 = ls_df[44, :scaled_load_shed][JJA_idx]
        lsJK = sum(ls_df[46:49, :load_shed])[JJA_idx]
        df[idx, :conditional_entropy] = estimate_ls_conditional_entropy(ls44)
        num_hrs_zoneJK_shortage = sum(lsJK .> 0.0)
        if num_hrs_zoneJK_shortage > 0
            # proportion of hours with load shedding in zone J or K in which there is also load shedding at bus 44
            df[idx, :prop_zoneJK_with_zoneG_shortage] = sum( (ls44 .> 0.0) .&& (lsJK .> 0.0) )/num_hrs_zoneJK_shortage
        else
            df[idx, :prop_zoneJK_with_zoneG_shortage] = Nothing
        end
    end
end
CSV.write(joinpath(output_dir, "scenario_features.csv"), df)