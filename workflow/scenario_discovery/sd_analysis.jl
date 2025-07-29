# activate environment
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

# load dependencies and utility functions
using CSV
using DataFrames
using DelimitedFiles
using MLJ

import Random
import StatsBase: countmap, nquantile, sample

include(joinpath(@__DIR__, "..", "utils", "other.jl"))

Random.seed!(11) # for reproducibility

# main script
# -----------

k = 3 # number of clusters

df = DataFrame(CSV.File(joinpath(@__DIR__, "..", "..", "output", "scenario_features.csv")))
labels = readdlm(joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "clustering", "cluster_labels_$(k).txt"), Int)
insertcols!(df, :label => vec(labels))

# transform features by discretizing
temp_edges = collect(1.5:1.5:7.5)
insertcols!(df, :disc_temp => map(x -> findfirst(x .<= temp_edges), df.avg_temp_anomaly))

wind_edges = collect(-3e3:1.5e3:4.5e3)
insertcols!(df, :disc_wind => map(x -> findfirst(x .<= wind_edges), df.avg_wind_anomaly))

solar_edges = collect(-3e3:1.5e3:6.0e3)
insertcols!(df, :disc_solar => map(x -> findfirst(x .<= solar_edges), df.avg_solar_anomaly))

hydro_edges = collect(-2.5e5:5.0e4:1.5e5)
insertcols!(df, :disc_hydro => map(x -> findfirst(x .<= hydro_edges), df.avg_hydro_anomaly))

# get right types
coerce!(df, :label       => Multiclass{k},
            :disc_temp   => OrderedFactor,
            :disc_wind   => OrderedFactor,
            :disc_solar  => OrderedFactor,
            :disc_hydro  => OrderedFactor,
            :batt_cap_sf => Continuous
)

# split into target and features
# we drop avg_load_anomaly because almost perfectly correlated with avg_temp_anomaly
y, X = unpack(df, ==(:label), in([:disc_temp, :disc_wind, :disc_solar, :disc_hydro, :batt_cap_sf]))

# hyperparameter tuning (:min_purity_increase, :merge_purity_threshold) and fit
Tree = @load DecisionTreeClassifier pkg=DecisionTree verbosity=0
tree = Tree(max_depth=3)
r1 = range(tree, :min_purity_increase, lower=0.001, upper=1.0, scale=:log)
r2 = range(tree, :merge_purity_threshold, lower=0.1, upper=1.0, scale=:linear)
self_tuning_tree = TunedModel(model=tree, 
                              resampling=StratifiedCV(nfolds=6), 
                              tuning=Grid(resolution=10), 
                              range=[r1, r2], 
                              measure=[brier_loss, cross_entropy, balanced_accuracy], 
                              repeats=3
)                              
mach = machine(self_tuning_tree, X, y)
fit!(mach)
rpt = report(mach)
display(rpt.best_history_entry.evaluation) # evaluation result for best hyperparameters
display(fitted_params(mach).best_fitted_params.tree)

# # fit basic decision tree (default hyperparameters)
# Tree = @load DecisionTreeClassifier pkg=DecisionTree verbosity=0
# tree = Tree(max_depth=3)
# mach = machine(tree, X, y)
# fit!(mach)

# examine results
println("--- CLUSTERING RESULTS ---")
display(confmat(mode.(predict(mach, X)), y)) # confusion matrix
bal_acc = balanced_accuracy(mode.(predict(mach, X)), y)
println("balanced accuracy: $(bal_acc)")

# save feature importances
CSV.write(joinpath(@__DIR__, "..", "..", "output", "scenario_discovery", "tree_feature_importance.csv"), stack(DataFrame(feature_importances(mach))))