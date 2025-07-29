# utility functions to set up and run independence tests for each source type (battery, generators, interface flows)

using Associations
using TimeseriesSurrogates: Surrogate
import StatsBase: Histogram, fit
import LinearAlgebra: normalize

include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))
include(joinpath(@__DIR__, "..", "utils", "independence_tests.jl"))

@kwdef struct IndependenceTestConfig
    k_embed::Int = 50    # number of NNs for target embedding optimization
    k_est::Int = 20      # number of NNs for kNN TE estimator
    bin_width = 0.2      # bin width for discretized TE estimator
    nshuffles::Int = 100 # number of shuffles
    w::Int = 23          # Theiler window
end

# utility function to retrieve pre-optimized embedding parameters
function get_embed_params(filename, scenario, yr, id, k)
    embed_params_df = DataFrame(CSV.File(joinpath(@__DIR__, "embedding_params", filename)))
    i = findall((embed_params_df.scenario .== scenario) .& 
                (embed_params_df.yr .== yr) .& 
                (embed_params_df.id .== id) .& 
                (embed_params_df.k .== k)
    )
    length(i) == 1 || throw(ArgumentError("Multiple entries for optimal embedding parameters."))

    return embed_params_df[i[1], :dim], embed_params_df[i[1], :delay]
end

# utility function to check if source time-series is (near-)constant so we can return nonsense result
function check_const_source(s, edges)
    h = fit(Histogram, s, edges)
    norm_h = normalize(h, mode=:probability)

    if maximum(norm_h.weights) > 0.97
        return true
    end

    return false
end

function batt_independence(load_bus::Int, batt::Int, scenario::Int, yr::Int, cfg::IndependenceTestConfig, surro_method::Surrogate)
    # load battery data
    src_df = scale_batt_results(scenario, yr)
    s_src = src_df[batt, :scaled_state][JJA_idx]

    if check_const_source(s_src, -0.01:0.01:1.01)
        return SurrogateAssociationTestResult(2, 0.0, zeros(cfg.nshuffles), 100.0, cfg.nshuffles)
    end
    
    # load load shed data
    ls_df = scale_load_shed_results(scenario, yr)
    s_trg = ls_df[load_bus, :scaled_load_shed][JJA_idx]

    # "optimize" target embedding - just reading here since we pre-optimized
    dT, τT = get_embed_params("ls_embed_params.csv", scenario, yr, load_bus, cfg.k_embed)

    # "optimize" source embedding - just reading here since we pre-optimized
    dS, τS = get_embed_params("batt_embed_params.csv", scenario, yr, batt, cfg.k_embed)

    # set up estimator
    emb = EmbeddingTE(dT=dT, τT=τT, dS=dS, τS=τS)
    meas = TEShannon(base=2, embedding=emb)
    est = CMIDecomposition(meas, FPVP(k=cfg.k_est, w=cfg.w))

    # set up independence test
    if (surro_method == MultidimRandomShuffle()) & (dS == 1)
        test = SurrogateAssociationTest(est, surrogate=RandomShuffle(), nshuffles=cfg.nshuffles)
    else
        test = SurrogateAssociationTest(est, surrogate=surro_method, nshuffles=cfg.nshuffles)
    end

    return independence(test, s_src.+1e-15randn(length(JJA_idx)), s_trg)
end

function flow_independence(load_bus::Int, flow::Int, scenario::Int, yr::Int, cfg::IndependenceTestConfig, surro_method::Surrogate)
    # load flow data
    src_df = scale_flow_results(scenario, yr)
    s_src = src_df[flow, :scaled_flow][JJA_idx]
    s_src = s_src ./ 2 # scale to have theoretical range ≤ 1

    if check_const_source(s_src, -0.51:0.01:0.51)
        return SurrogateAssociationTestResult(2, 0.0, zeros(cfg.nshuffles), 100.0, cfg.nshuffles)
    end
    
    # load load shed data
    ls_df = scale_load_shed_results(scenario, yr)
    s_trg = ls_df[load_bus, :scaled_load_shed][JJA_idx]
    
    # "optimize" target embedding - just reading here since we pre-optimized
    dT, τT = get_embed_params("ls_embed_params.csv", scenario, yr, load_bus, cfg.k_embed)

    # "optimize" source embedding - just reading here since we pre-optimized
    dS, τS = get_embed_params("flow_embed_params.csv", scenario, yr, flow, cfg.k_embed)

    emb = EmbeddingTE(dT=dT, τT=τT, dS=dS, τS=τS)
    meas = TEShannon(base=2, embedding=emb)

    # set up estimator
    if flow .>= 95 
        # dc lines more discrete-valued so using discretization + binning estimator
        precise = true
        disc = CodifyVariables(ValueBinning(RectangularBinning(cfg.bin_width, precise)))
        est = EntropyDecomposition(meas, PlugIn(Shannon()), disc)
    else
        s_src .+= 1e-15randn(length(JJA_idx))
        est = CMIDecomposition(meas, FPVP(k=cfg.k_est, w=cfg.w))
    end
    
    # set up independence test
    if (surro_method == MultidimRandomShuffle()) & (dS == 1)
        test = SurrogateAssociationTest(est, surrogate=RandomShuffle(), nshuffles=cfg.nshuffles)
    else
        test = SurrogateAssociationTest(est, surrogate=surro_method, nshuffles=cfg.nshuffles)
    end

    return independence(test, s_src, s_trg)
end

function if_flow_independence(load_bus::Int, if_flow::Int, scenario::Int, yr::Int, cfg::IndependenceTestConfig, surro_method::Surrogate)
    # load if_flow data
    src_df = scale_if_flow_results(scenario, yr)
    s_src = src_df[if_flow, :scaled_if_flow][JJA_idx]
    s_src = (s_src ./ 2) # scale to have theoretical range ≤ 1
    
    if check_const_source(s_src, -0.51:0.01:0.51)
        return SurrogateAssociationTestResult(2, 0.0, zeros(cfg.nshuffles), 100.0, cfg.nshuffles)
    end

    # load load shed data
    ls_df = scale_load_shed_results(scenario, yr)
    s_trg = ls_df[load_bus, :scaled_load_shed][JJA_idx]

    # "optimize" target embedding - just reading here since we pre-optimized
    dT, τT = get_embed_params("ls_embed_params.csv", scenario, yr, load_bus, cfg.k_embed)

    # "optimize" source embedding - just reading here since we pre-optimized
    dS, τS = get_embed_params("if_flow_embed_params.csv", scenario, yr, if_flow, cfg.k_embed)

    # set up estimator
    emb = EmbeddingTE(dT=dT, τT=τT, dS=dS, τS=τS)
    meas = TEShannon(base=2, embedding=emb)
    est = CMIDecomposition(meas, FPVP(k=cfg.k_est, w=cfg.w))
    
    # set up independence test
    if (surro_method == MultidimRandomShuffle()) & (dS == 1)
        test = SurrogateAssociationTest(est, surrogate=RandomShuffle(), nshuffles=cfg.nshuffles)
    else
        test = SurrogateAssociationTest(est, surrogate=surro_method, nshuffles=cfg.nshuffles)
    end

    return independence(test, s_src.+1e-15randn(length(JJA_idx)), s_trg)
end

function gen_independence(load_bus::Int, gen::Int, scenario::Int, yr::Int, cfg::IndependenceTestConfig, surro_method::Surrogate)
    # load gen data
    src_df = scale_gen_results(scenario, yr)
    s_src = src_df[gen, :scaled_gen][JJA_idx]

    # load load shed data
    ls_df = scale_load_shed_results(scenario, yr)
    s_trg = ls_df[load_bus, :scaled_load_shed][JJA_idx]

    # "optimize" target embedding - just reading here since we pre-optimized
    dT, τT = get_embed_params("ls_embed_params.csv", scenario, yr, load_bus, cfg.k_embed)

    # "optimize" source embedding - just reading here since we pre-optimized
    if (src_df[gen, :type] == "Import") & (gen ∉ [11, 16])
        # imported generation (except gens 11 and 16 - these gens have more interesting dynamics)
        dS = 1
        τS = -1
    else
        dS, τS = get_embed_params("gen_embed_params.csv", scenario, yr, gen, cfg.k_embed)
    end

    # set up estimator
    emb = EmbeddingTE(dT=dT, τT=τT, dS=dS, τS=τS)
    meas = TEShannon(base=2, embedding=emb)
    precise = true
    disc = CodifyVariables(ValueBinning(RectangularBinning(cfg.bin_width, precise)))
    est = EntropyDecomposition(meas, PlugIn(Shannon()), disc)

    # set up independence test
    if (surro_method == MultidimRandomShuffle()) & (dS == 1)
        test = SurrogateAssociationTest(est, surrogate=RandomShuffle(), nshuffles=cfg.nshuffles)
    else
        test = SurrogateAssociationTest(est, surrogate=surro_method, nshuffles=cfg.nshuffles)
    end

    return independence(test, s_src, s_trg)
end

function curtail_independence(load_bus::Int, curtail::Int, scenario::Int, yr::Int, cfg::IndependenceTestConfig, surro_method::Surrogate)
    # load curtail data
    src_df = scale_curtail_results(scenario, yr)
    s_src = src_df[curtail, :scaled_curtail][JJA_idx]

    # load load shed data
    ls_df = scale_load_shed_results(scenario, yr)
    s_trg = ls_df[load_bus, :scaled_load_shed][JJA_idx]

    # "optimize" target embedding - just reading here since we pre-optimized
    dT, τT = get_embed_params("ls_embed_params.csv", scenario, yr, load_bus, cfg.k_embed)

    # "optimize" source embedding - just reading here since we pre-optimized
    dS, τS = get_embed_params("curtail_embed_params.csv", scenario, yr, curtail, cfg.k_embed)

    # set up estimator
    emb = EmbeddingTE(dT=dT, τT=τT, dS=dS, τS=τS)
    meas = TEShannon(base=2, embedding=emb)
    precise = true
    disc = CodifyVariables(ValueBinning(RectangularBinning(cfg.bin_width, precise)))
    est = EntropyDecomposition(meas, PlugIn(Shannon()), disc)

    # set up independence test
    if (surro_method == MultidimRandomShuffle()) & (dS == 1)
        test = SurrogateAssociationTest(est, surrogate=RandomShuffle(), nshuffles=cfg.nshuffles)
    else
        test = SurrogateAssociationTest(est, surrogate=surro_method, nshuffles=cfg.nshuffles)
    end

    return independence(test, s_src, s_trg)
end