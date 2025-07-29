# functions to optimize embedding parameters according to the Ragwitz criterion as described in Ragwitz and Kantz (2002)
# based on TRENTOOL implementation (https://github.com/trentool/TRENTOOL3/blob/ed2d45b744c91693632fb6bcf8a6301acb5ea2a2/private/TEragwitz.m#L29)

using DelayEmbeddings
using Distances
using Neighborhood
using Statistics: std

# struct to store parameters for embedding optimization
@kwdef struct EmbeddingOptions
    max_d::Int # max embedding dimension

    max_τ::Int # max embedding delay

    nbhd_type::String
    # nbhd_type: 
    # - "mass" uses NeighborNumber(nbhd_param)
    # - "range" uses WithinRange(nbhd_param)

    nbhd_param::Union{Int, Float64}
    # nbhd_param:
    # - if nbhd_type = "mass", number of NNs to use in prediction
    # - if nbhd_type = "range", distance within which to find NNs

    u::Int = 1 # prediction lag

    w::Int = 0 # Theiler window

    metric::Metric = Cityblock() # metric for finding NNs in state space
end

# predicts time series using local model
function find_local_pred(s::Vector, query_idx::Vector{Int}, d::Int, τ::Int, opt::EmbeddingOptions)
    # check nbhd_type is valid
    opt.nbhd_type ∈ ["mass", "range"] || throw(ArgumentError("nbhd_type must be either 'mass' or 'range'"))
    
    # check queries are in appropriate range
    maximum(query_idx) ≤ length(s) || throw(ArgumentError("elements of query_idx must be ≤ $(length(s))"))
    minimum(query_idx) ≥ 1+(opt.max_d-1)*opt.max_τ+opt.u || throw(ArgumentError("elements of query_idx must be ≥ $(1+(opt.max_d-1)*opt.max_τ+opt.u)"))
    n_queries = length(query_idx)
    
    # convert query idxs to predictor idxs
    pred_idx = query_idx .- ((opt.max_d-1)*opt.max_τ + opt.u)

    # embed time series
    noisy_s = s .+ 1e-15randn(length(s)) # add noise to avoid numerical error with NN search
    ssset = genembed(noisy_s[1+(opt.max_d-1)*opt.max_τ-(d-1)*τ : end-opt.u], collect(0:-τ:-(d-1)*τ))
    # effective length of s is length(ssset) = length(s) - (max_d - 1)*max_τ - u
    # so length(ssset) consistent across all (d,τ)-combos

    # preprocessing for nearest neighbor search
    tree = searchstructure(KDTree, ssset, opt.metric)
    theiler = Theiler(opt.w, pred_idx)
    
    if opt.nbhd_type == "mass"
        # find indices of nearest nbhd_param neighbors
        nbhd_idx = bulkisearch(tree, ssset[pred_idx], NeighborNumber(opt.nbhd_param), theiler)
        # nbhd_idx is vector of vectors

        # create local predictions
        local_pred = zeros(n_queries)
        for i = 1:n_queries
            local_pred[i] = sum(noisy_s[nbhd_idx[i] .+ (opt.u+(opt.max_d-1)*opt.max_τ)]) / opt.nbhd_param
        end
    else
        # find indices of all neighbors at distance ≤ nbhd_param 
        nbhd_idx = bulkisearch(tree, ssset[pred_idx], WithinRange(opt.nbhd_param), theiler)
        counts = length.(nbhd_idx)
        
        # make sure all neighborhoods are nonempty
        sum(counts .== 0) == 0 || throw(ArgumentError("range too small, found neighborhood with no neighbors"))

        # create local predictions
        local_pred = zeros(n_queries)
        for i = 1:n_queries
            local_pred[i] = sum(noisy_s[nbhd_idx[i] .+ (opt.u+(opt.max_d-1)*opt.max_τ)]) / counts[i]
        end
    end

    return local_pred
end

# calculates MRE for local predictions with embedding dimension d and delay time τ
function ragwitz_mre(s::Vector, query_idx::Vector{Int}, d::Int, τ::Int, opt::EmbeddingOptions)
    # create local predictions
    local_pred = find_local_pred(s, query_idx, d, τ, opt)

    # calculate MRE and normalize by std of series
    return ragwitz_mre(s, query_idx, local_pred)
end

# calculates MRE for local predictions with embedding dimension d and delay time τ
function ragwitz_mre(s::Vector, query_idx::Vector{Int}, local_pred::AbstractVector)
    length(query_idx) == length(local_pred) || throw(DimensionMismatch("length(query_idx) must match length(local_pred)"))

    return ( sum( (local_pred - s[query_idx]).^2 ) / length(query_idx) )/std(s[query_idx])
end

# calculates misclassification error for local predictions
function ragwitz_misclass_err(s::Vector, query_idx::Vector{Int}, d::Int, τ::Int, opt::EmbeddingOptions)
    # create local predictions
    local_pred = find_local_pred(s, query_idx, d, τ, opt)

    return ragwitz_misclass_err(s, query_idx, local_pred)
end

# calculates misclassification error for local predictions
function ragwitz_misclass_err(s::Vector, query_idx::Vector{Int}, local_pred::Vector)
    length(query_idx) == length(local_pred) || throw(DimensionMismatch("length(query_idx) ($(length(query_idx))) must match length(local_pred) ($(length(local_pred)))"))
    
    s01 = (s[query_idx] .> 1e-8)
    pred01 = (local_pred .> 0.05)

    # calculate fraction of timesteps where prediction doesn't match series
    return 1 - sum( s01 .== pred01 ) / length(query_idx)
end

# returns optimal embedding parameters
function optimize_embedding(s::Vector, opt::EmbeddingOptions; query_idx::Vector{Int}=collect(1+(opt.max_d-1)*opt.max_τ+opt.u : length(s)), n_itr=10)
    err = zeros(opt.max_d, opt.max_τ)

    for d in 1:opt.max_d, τ in 1:opt.max_τ
        for i in n_itr
            local_pred = find_local_pred(s, query_idx, d, τ, opt)
            err[d, τ] += ragwitz_mre(s, query_idx, local_pred)
        end
    end

    # find d and τ corresponding to minimum error
    d = findmin(err)[2][1]
    if d == 1
        return d, -1 # if d = 1, embedded "vector" is just value at previous timestep
    else
        return d, -findmin(err)[2][2]
    end
end