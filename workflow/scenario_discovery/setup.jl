# utility functions to set up and run independence tests for IF flows

using Associations

include(joinpath(@__DIR__, "..", "utils", "embedding.jl"))
include(joinpath(@__DIR__, "..", "utils", "scaling.jl"))
include(joinpath(@__DIR__, "..", "utils", "independence_tests.jl"))

@kwdef struct IndependenceTestConfig
    k_embed::Int = 50        # number of NNs for target embedding optimization
    k_est::Int = 20          # number of NNs for kNN TE estimator
    bins::AbstractRange      # bins for discretized TE estimator
    nshuffles::Int = 100     # number of shuffles
    w::Int = 23              # Theiler window
end

function if_flow_independence(s_trg, if_flow::Int, scenario::Int, yr::Int, cfg::IndependenceTestConfig, dT::Int, τT::Int)
    # load IF flow data
    src_df = scale_if_flow_results(scenario, yr)
    s_src = src_df[if_flow, :scaled_if_flow][JJA_idx] ./ 2 # scale to have theoretical range ≤ 1

    # optimize source embedding
    opt = EmbeddingOptions(max_d=8, max_τ=24, nbhd_type="mass", nbhd_param=cfg.k_embed, w=23, metric=Cityblock())
    dS, τS = optimize_embedding(s_src, opt; n_itr=20)
    
    # set up independence test
    emb = EmbeddingTE(dT=dT, τT=τT, dS=dS, τS=τS) # target embedding was passed in as argument
    meas = TEShannon(base=2, embedding=emb)
    if if_flow ∈ [6, 13]
        # shift s_src to fall in [0, 1] for more straightforward discretization 
        s_src .+= 0.5
        s_src[s_src .<  0.0] .= 0.0

        # set up discretization and discrete estimator      
        precise = true
        disc = CodifyVariables(ValueBinning(FixedRectangularBinning(cfg.bins, 1, precise)))
        est = EntropyDecomposition(meas, PlugIn(Shannon()), disc)
        test = SurrogateAssociationTest(est, surrogate=WLS(IAAFT(), rescale=true), nshuffles=cfg.nshuffles)

        return independence(test, s_src, s_trg)
    else
        # set up kNN estimator
        est = CMIDecomposition(meas, FPVP(k=cfg.k_est, w=cfg.w))
        test = SurrogateAssociationTest(est, surrogate=WLS(IAAFT(), rescale=true), nshuffles=cfg.nshuffles)

        # add a bit of noise to time series to help break NN ties
        return independence(test, s_src+1e-15randn(length(JJA_idx)), s_trg+1e-15randn(length(JJA_idx)))
    end
end