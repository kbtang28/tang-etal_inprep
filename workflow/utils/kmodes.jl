# implementation of K-modes clustering algorithm as described in Huang (1998) and Ng et al. (2007)

using Printf
import Random
import StatsBase

mutable struct KmodesResult{M <: AbstractMatrix{Int}, D <: Real}
    modes::M                    # cluster modes (d x k matrix)
    assignments::Vector{Int}    # assignments (n)
    counts::Vector{Int}         # number of points assigned to each cluster (k)
    total_cost::D               # objective
    iterations::Int             # number of elapsed iterations
    converged::Bool             # whether algorithm converged
end

function matching_dissim(a, b)
    return sum(a .!= b)
end

function ng_dissim(a::AbstractMatrix{Int}, b::Vector{Int}, W::AbstractMatrix{Bool}, data::AbstractMatrix{Int})
    # a    = matrix of k modes (d x k)
    # b    = query point (d x 1)
    # W    = membership matrix (k x n)
    # data = data matrix (d x n)

    function calc_cljr(b, data, meml::Vector{Bool}, j::Int)
        return sum(data[j, meml] .== b[j])
    end

    function calc_φ(b, data, meml::Vector{Bool}, j::Int)
        cl = sum(meml) # size of l-th cluster
        return cl > 0 ? 1.0 - (calc_cljr(b, data, meml, j) / cl) : 0.0
    end

    size(a, 2) == size(W, 1) || throw(DimensionMismatch("Got $(size(a,2)) modes in a but $(size(W, 1)) clusters in W"))
    size(W, 2) == size(data, 2) || throw(DimensionMismatch("Got cluster membership matrix for $(size(W, 2)) points but data matrix has $(size(data, 2)) points"))

    return [sum([b[j] == r ? calc_φ(b, data, W[i, :], j) : 1.0 for (j, r) in enumerate(vals_a)]) for (i, vals_a) in enumerate(eachcol(a))]
end

function kmodes(data::AbstractMatrix{Int}, k::Int; seed::Union{Nothing, AbstractMatrix{Int}}=nothing, max_iter=200, verbose=false)
    n = size(data, 2)
    
    # initialize modes
    if seed == nothing
        unique_ids = unique(i -> data[:, i], 1:n)
        modes = data[:, StatsBase.sample(unique_ids, k, replace=false)]
    else
        (size(seed, 1) == size(data, 1)) && (size(seed, 2) == k) || throw(ArgumentError("Seeds of wrong dimension: got $(size(seed, 1)) x $(size(seed, 2)) matrix, need $(size(data, 1)) x $(k) matrix"))
        modes = seeds
    end

    # assignment matrix
    W = zeros(Bool, k, n) 

    # initial assignments to clusters
    @inbounds for j in 1:n
        a = argmin(ng_dissim(modes, data[:, j], W, data)) # note W updated with each point processed!
        W[a, j] = true
    end

    # initial mode update
    _kmodes_find_modes!(modes, data, W)

    # main loop
    t = 0
    converged = false
    obj = _kmodes_get_obj(modes, data, W)

    if verbose
        @printf("%7s %18s %18s %18s\n", "iter", "obj", "obj-change", "n_changes")
        println("-----------------------------------------------------------------------")
        @printf("%7d %18.6e\n", t, obj)
    end

    while !converged && (t < max_iter)
        t += 1

        # one iteration
        n_changes = _kmodes_iter!(modes, data, W)
        
        pobj = obj
        obj = _kmodes_get_obj(modes, data, W)

        # check convergence
        converged = ((n_changes == 0) || (obj >= pobj))

        if verbose
            @printf("%7d %18.6e %18.6e %18d\n", t, obj, obj - pobj, n_changes)
        end
    end

    # final assignment to clusters (no updating modes)
    @inbounds for j in 1:n
        a = argmin(ng_dissim(modes, data[:, j], W, data)) # note W updated with each point processed!
        pa = findall(view(W, :, j))[1]
        if pa == a
            continue
        else
            W[pa, j] = false
            W[a,  j] = true
        end
    end
    obj = _kmodes_get_obj(modes, data, W)

    if verbose
        if converged
            println("K-modes converged with $(t) iterations (objv = $(obj))")
        else
            println("K-modes terminated without convergence after $(t) iterations (objv = $(obj))")
        end
    end

    # return results
    assignments = map(w -> findall(w)[1], eachcol(W))
    counts = sum(W, dims=2)
    return KmodesResult(modes, assignments, vec(counts), obj, t, converged)
end

function _kmodes_iter!(modes::AbstractMatrix{Int}, data::AbstractMatrix{Int}, W::AbstractMatrix{Bool})
    n = size(data, 2)
    n_changes = 0

    @inbounds for j in 1:n
        a = argmin(ng_dissim(modes, data[:, j], W, data))
        pa = findall(view(W, :, j))[1]
        if pa == a
            # pt already in the right place
            continue
        else 
            n_changes += 1

            # reassign
            W[pa, j] = false
            W[a,  j] = true

            # update mode of new cluster
            modes[:, a] .= mapslices(x -> minimum(StatsBase.modes(x)), data[:, W[a, :]], dims=2)

            # update mode of old cluster
            if sum(W[pa, :]) == 0
                # if old cluster empty, reinitialize with random point from largest cluster
                l = argmax(sum(W, dims=2))
                modes_set = Set(eachcol(modes))
                for i in Random.shuffle(findall(W[l, :]))
                    if data[:, i] ∉ modes_set
                        modes[:, pa] .= data[:, id]
                        break
                    end
                end
            else
                modes[:, pa] .= mapslices(x -> minimum(StatsBase.modes(x)), data[:, W[pa, :]], dims=2)
            end
        end
    end
    
    @assert sum(W) == n "Uh oh... one or more data points have been improperly assigned!"
    return n_changes
end

function _kmodes_find_modes!(modes::AbstractMatrix{Int}, data::AbstractMatrix{Int}, W::AbstractMatrix{Bool})
    d, nD = size(data)
    k, nW = size(W)

    nD == nW || throw(DimensionMismatch("Got cluster membership matrix for $(nW) points but data matrix has $(nD) points"))

    for l in 1:k
        meml = W[l, :]
        if sum(meml) == 0
            # if cluster is empty, choose new mode from data at random
            modes_set = Set(eachcol(modes))
            for i in Random.shuffle(unique(i -> data[:, i], 1:nD))
                if data[:, i] ∉ modes_set
                    modes[:, l] .= data[:, i]
                    break
                end
            end
        else
            # each feature is set to mode for current cluster
            modes[:, l] .= mapslices(x -> minimum(StatsBase.modes(x)), data[:, meml], dims=2)
        end
    end
end

function _kmodes_get_obj(modes::AbstractMatrix{Int}, data::AbstractMatrix{Int}, W::AbstractMatrix{Bool})
    n = size(data, 2)
    obj = 0.0

    @inbounds for j in 1:n
        obj += ng_dissim(modes, data[:, j], W, data)[W[:, j]][1]
    end

    return obj
end
