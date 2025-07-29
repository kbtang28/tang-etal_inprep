_your zenodo badge here_

# Tang-etal_inprep

**Identification of pressure points in decarbonized power systems using transfer entropy**

Katerina Tang<sup>1\*</sup>, M. Vivienne Liu<sup>2</sup>, C. Lindsay Anderson<sup>1,3,</sup>, and Vivek Srikrishnan<sup>3</sup>

<sup>1</sup> Center for Applied Mathematics, Cornell University, Ithaca, NY

<sup>2</sup> Systems Engineering, Cornell University, Ithaca, NY

<sup>3</sup> Department of Biological & Environmental Engineering, Cornell University, Ithaca, NY

\* corresponding author:  kbt28@cornell.edu

## Abstract

Renewable energy integration and end-use electrification increase the weather-dependence of decarbonized energy systems. 
Identifying critical infrastructure that constrains the power grid's ability to meet electricity demand under weather-induced shocks and stressors is essential for understanding vulnerabilities and guiding adaptation. 
We use transfer entropy to identify predictive pressure points: grid components whose utilization patterns provide early signals of downstream power shortages.
Applied to simulations of New York State's proposed future zero-carbon grid under a range of meteorological and technological scenarios, our method shows that these pressure points often arise from complex, system-wide interactions between generation, transmission, and demand.
While transfer entropy does not support strong conclusions about causality, the identified pressure points align with known vulnerabilities and offer insight into failure pathways.
Furthermore, these pressure points are not easily predicted by high-level scenario features alone, underscoring the need for holistic and adaptive approaches to reliability planning in low-carbon power systems.

## Journal reference
Link to preprint: incoming!

## Data reference

### Input data & contributing modeling software

We use simulations from [ACORN](https://github.com/AndersonEnergyLab-Cornell/ny-clcpa2050/tree/main) in our study. For convenience, relevant model input and simulation results for the three base scenarios analyzed in the first part of our results are provided in the `ACORN` directory. Reproducing the scenario discovery portion of our analysis requires model input and simulation results from the full scenario ensemble described in [Liu et al. (2023)](https://arxiv.org/abs/2307.15079) and documented in the above repository.

### Output data

Output from our analysis can be found in the `output` directory.

- `output/base_scenarios` has summary results of the independence tests for all three base scenarios.
- `output/scenario_discovery` has summary results of the independence tests for the entire scenario ensemble and results from the cluster analysis.
- `scenario_features.csv` has relevant meteorological and technological features for each scenario-year pair. This data is used in the scenario discovery analysis and to create some figures.

## Dependencies

This code is based on Julia 1.10.0. Relecant dependencies are in the `Project.toml` and `Manifest.toml` files. (The `Manifest.toml` file specifies the particular versions; this file should be kept as-is for perfect reproducibility but may need to be deleted and rebuilt with `Pkg.instantiate()` for different Julia versions.)

## Reproduction

After cloning the repository, install the necessary packages:
```julia
import Pkg
Pkg.activate(".") # from cloned root directory
Pkg.instantiate()
```

To summarize the meteorological and technological features for each scenario-year pair, run `julia workflow/scenario_features.jl`.

### Base scenarios
Scripts to re-run the analysis for the three base scenarios are in the `workflow/base_scenarios/` directory.
| Script name | Description & notes | How to run (from root directory)|
| --- | --- | --- |
| `compare_embedding_opts.jl` | Script to pre-optimize embedding parameters for TE estimates. Saves results and diagnostic figures in `workflow/base_scenarios/embedding_params/`. | `sbatch workflow/base_scenarios/compare_embedding_opts.sh` |
| `generate_wls_surrogates.jl` | Script to generate WLS surrogates. Saves results in `workflow/base_scenarios/precomputed_wls_surrogates/`. | `sbatch workflow/base_scenarios/generate_wls_surrogates.sh` |
| `test_independence.jl` | Script to perform pairwise TE analysis. Saves results in `output/base_scenarios/`. Setup for independence tests are in `setup.jl`. | `sbatch workflow/base_scenarios/test_independence.sh` |


### Scenario discovery
Scripts to re-run the scenario discovery analysis are in the `workflow/scenario_discovery/` directory.
| Script name | Description & notes | How to run (from root directory)|
| --- | --- | --- |
| `test_independence.jl` | Script to perform pairwise TE analysis. Saves results in `output/scenario_discovery/independence/`. Setup for independence tests are in `setup.jl`. | `sbatch workflow/scenario_discovery/test_independence.sh` |
| `cluster_analysis.jl` | Script to cluster results from independence tests. Saves results in `output/scenario_discovery/clustering/`. | `sbatch workflow/scenario_discovery/cluster_analysis.sh` |
| `sd_analysis.jl` | Script to fit classification tree to predict cluster membership. | `julia workflow/scenario_discovery/sd_analysis.jl` |

## Reproduce paper figures
The scripts listed below reproduce the named figures. Output files correspond to the files in `figures/`.

| Figure(s) | Script | Output File |
| --- | --- | --- |
| Fig. 2  | `workflow/base_scenarios/plot_base_scenarios.jl` | `scenario_features.png` |
| Fig. 3 |`workflow/base_scenarios/map_pressure_points.ipynb` | `map_pressure_points.png` |
| Fig. 4 | `workflow/base_scenarios/plot_base_scenarios.jl` | `rewnew_ratios.png` |
| Fig. 5 | `workflow/scenario_discovery/map_simplified_IFs.ipynb` | `map_3_clusters.png` |
| Fig. 6 | `workflow/scenario_discovery/plot_temp_and_solar_by_cluster.jl` | `clusters_solar_temp_boxplots.png` |
| Fig. S1 | `workflow/basic_te_example.jl` | `basic_te_example.png` |
| Fig. S2 | `workflow/scenario_discovery/map_zonal_caps.ipynb` | `map_zonal_caps.png` |
| Fig. S3 | `workflow/plot_load_and_hydro_vs_temp.jl` | `load_and_hydro_vs_temp.png`|
| Fig. S4 | `workflow/plot_ls_by_bus.jl` | `ls_hrs_and_prop_by_bus.png` |
| Fig. S5 | `workflow/base_scenarios/plot_base_scenarios.jl` | `curtailment.png` |
| Fig. S6 | `workflow/base_scenarios/plot_well_behaved_scenario.jl` | `s140_GHif_zoneFrenewables.png` |
| Fig. S7 | `workflow/base_scenarios/plot_limited_resource_scenario.jl` | `s69_ABif_BCif_IJif_utilization.png`|
| Fig. S8 | `workflow/base_scenarios/plot_limited_resource_scenario.jl` | `s69_zoneJwind_IJif_scatter.png`|
| Fig. S9 | `workflow/base_scenarios/plot_extreme_temp_scenario.jl` | `s290_ABif_CEif_GHif_HIif_IKif_utilization.png` |
| Fig. S10 | `workflow/base_scenarios/plot_extreme_temp_scenario.jl` | `s290_HIif_IKif_utilization.png` |
| Fig. S11 | `workflow/scenario_discovery/plot_example_IF_dynamics.jl` | `clus1_ABif_BCif_PJMif_utilization.png` |
| Fig. S12 | `workflow/scenario_discovery/plot_IKif_dynamics.jl` | `clus1_IKif_util.png` |
| Fig. S13 | `workflow/plot_EGif_examples.jl` | `EGif_examples.png` |
| Fig. S14 | `workflow/scenario_discovery/cluster_analysis.jl` | `cluster_obj_vs_k.png` |
| Fig. S15 | `workflow/scenario_discovery/map_simplified_IFs.ipynb` | `map_3_clusters.png` |
| Fig. S17 | `workflow/plot_embed_err_examples.jl` | `ragwitz_criterion.png` |
| Fig. S18 | `workflow/plot_wls_surrogates.sh` | `wls_surro_examples.png` |