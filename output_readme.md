# readme for output files

the `output` folder contains files:

```
output/
├── g1-g4/
 ├── g4_hwmatches.csv
 ├── h1_byparamid_filtered.csv
 ├── h1_byrep_filtered.csv
 ├── h4_byparamid_filtered.csv
 ├── h4_byrep_filtered.csv
 └── h4_type1_fromwr.csv
├── popnets/
 ├── topgraphs_1graph_byirep-june27.csv
 ├── topgraphs_1graph_byparamset-june27.csv
 ├── topgraphs_3thresh_byirep-june27.csv
 ├── topgraphs_3thresh_byparamset-june27.csv
 ├── topgraphs_5graphs_byirep-june27.csv
 ├── topgraphs_5graphs_byparamset-june27.csv
 ├── topgraphs_byparamset_filtered-june27.csv
 └── topgraphs_byrep_filtered-june27.csv
├── citations_clean.csv
├── nature_science.csv
└── timingdata.csv
```

## g1-g4

### g4_hwmatches.csv

This file contains the graph search "matches" for the `g4` graph,
where they did not match based on `hash` (`truenet == FALSE`),
but do have hardwired cluster distance of 0 (`net_hwdist == 0`).
Generated in [`summarize_graphs.qmd`](../scripts/summarize_graphs.qmd).

**g4 hw matches data column structure**:

- `newick`: newick string of estimated network
- `score`: log-likelihood score for a graph as given by admixtools
- `hash`: unique value given to a given graph topology by admixtools. "for two graphs to be identical, the topology and leaf nodes have to match, but internal node names do not matter." used to check for isomorphism.
- `truenet`: boolean true/false if graph is a match to the true net based on `hash`
- `net_hwdist`: hardwired cluster distance from true net (0 if net is true net)
- `tree_hwdist`: minimum hardwired cluster distance from displayed trees of inferred graph to major tree of true graph
- `scorediff`: `truenet`- "estimated score"'s score (0 if net is true net)
        this is the loglikelihood of the estimated graph minus the loglikelihood
        of the true topology. higher is better!
- `h_search`: number of max number of hybridization events searched for (`max_admix` argument to `find_graphs`)
- `h_est`: number of hybridization events estimated
- `ndisplayed`: number of displayed trees from inferred graph, with gamma threshold of 0.0,
        that match near-major trees from true network
- `irep`: replicate number
- `irun`: run number
- `lindist`: value for standard deviation on the log scale of log-normal distribution with mean 1, used to pull r_l values.
- `nind`: number of individuals per population
- `nbiallelic`: number of biallelic sites simulated
- `mutrate`: mutation rate simulated
- `net`: which of the 18 fixed starting networks data was simulated under
- `paramid`: string combination of all the possible parameters

### h1 & h4 by rep

Generated in [`summarize_graphs.qmd`](../scripts/summarize_graphs.qmd).

.csv file containing results for graph search of g1 (`h1_byrep_filtered.csv`) or g4 (`h4_byrep_filtered.csv`) networks,
as summarized by replicate.

**h1/h4 by replicate data column structure**:

- `irep`: replicate number
- `paramid`: string combination of all the possible parameters
- `h_search`: number of max number of hybridization events searched for (`max_admix` argument to `find_graphs`)
- `lindist`: value for standard deviation on the log scale of log-normal distribution with mean 1, used to pull r_l values
- `nind`: number of individuals per population
- `net`: which of the 18 fixed starting networks data was simulated under
- `mutrate`: mutation rate simulated
- `nbiallelic`: number of biallelic sites simulated
- `ngraphs`: sum of graphs retained across the replicate
- `min_nethw`: minimum hardwired cluster distance from inferred graphs to true graph
- `atleastonehashmatch`: boolean value whether at least 1 graph that matches the `hash` value of true graph (isomorphic match)
- `min_treehw`: minimum (across replicate) hardwired cluster distance from displayed trees of inferred graph to major tree of true graph
- `max_ndisplayed`: maximum number of displayed trees from inferred graph, with gamma threshold of 0.0, that match either of 2 near-major trees from true network (gamma > 0.45). values should not exceed 2. (0 = bad, 2 = best)
- `mean_ndisplayed`: mean number of displayed trees from inferred graph, with gamma threshold of 0.0, that match either of 2 near-major trees from true network (gamma > 0.45). values should not exceed 2. (0 = bad, 2 = best)
- `relativell`: maximum value of `scorediff` (= score of truenet - score of estimated net), the best score of all inferred graphs.

### h1 & h4 by paramid

Generated in [`summarize_graphs.qmd`](../scripts/summarize_graphs.qmd).

.csv file containing results for graph search of g1 (`h1_byparamid_filtered.csv`) or g4 (`h4_byparamid_filtered.csv`) networks, as summarized by parameter set.

**h1/h4 by parameter set data column structure**:

- `irep`: replicate number
- `paramid`: string combination of all the possible parameters
- `h_search`: number of max number of hybridization events searched for (`max_admix` argument to `find_graphs`)
- `lindist`: value for standard deviation on the log scale of log-normal distribution with mean 1, used to pull r_l values
- `nind`: number of individuals per population
- `net`: which of the 18 fixed starting networks data was simulated under
- `mutrate`: mutation rate simulated
- `nbiallelic`: number of biallelic sites simulated
- `prop_min_nethw_eq_0`: proportion of replicates that have at least 1 graph with hardwired cluster distance of 0 to true network
- `prop_min_hashmatch`: proportion of replicates that have at least 1 graph that matches the `hash` value of true graph (isomorphic match)
- `prop_min_treehw_eq_0`: proportion of replicates that have at least 1 graph with hardwired cluster distance of 0 to major tree of true graph
- `ngraphs`: sum of graphs retained across the parameter set
- `prop_max_ndisp_eq_2`: proportion of replicates that have at least 1 graph that displays both near major trees of the true network.

### h1_type1_fromwr.csv and h4_type1_fromwr.csv

Summary of results when evaluating type-1 error from using threshold that graphs must have
a worst residual no greater than 3.
Generated in [`summarize_graphs.qmd`](../scripts/summarize_graphs.qmd) from data produced in
`calcwr_g1.jl` (results in `h1_type1_fromwr.csv`) and `calcwr_g4.jl`
(results in `h4_type1_fromwr.csv`).

**type 1 error from worst residuals data column structure in both files**:

- `paramid`: string combination of all the possible parameters
- `lindist`: value for standard deviation on the log scale of log-normal distribution with mean 1, used to pull r_l values
- `net`: which of the 18 fixed starting networks data was simulated under
- `mutrate`: mutation rate simulated
- `nbiallelic`: number of biallelic sites simulated
- `type1`: type 1 error rate (1 - `prop_whsearch1` or `prop_whsearch4` - so proportion of replicates who have zero graphs remaining after applying threshold of retaining graphs with worst residual values no greater than 3)

**type 1 error from worst residuals data column structure in `h1_type1_fromwr.csv`**:

- `prop_whsearch1`: proportion of replicates who have a non-zero amount of graphs after applying threshold of retaining graphs with worst residual absolute values no greater than 3 when `h_search == 1`.
- `median_ngraphs_h1_before`: median number of graphs across across 100 replicates retained when `h_search == 1` before applying threshold of retaining graphs with worst residual absolute values no greater than 3.
- `median_ngraphs_h1_after`: median number of graphs across across 100 replicates retained when `h_search == 1` after applying threshold of retaining graphs with worst residual absolute values no greater than 3.
- `median_ngraphs_h2_before`: same as `median_ngraphs_h1_before`, but when `h_search == 2`
- `median_ngraphs_h2_after`: same as `median_ngraphs_h1_after`, but when `h_search == 2`

**type 1 error from worst residuals data column structure in `h4_type1_fromwr.csv`**:

- `prop_whsearch4`: proportion of replicates who have a non-zero amount of graphs after applying threshold of retaining graphs with worst residual absolute values no greater than 3 when `h_search == 4`.
- `median_ngraphs_h4_before`: median number of graphs across 100 replicates retained when `h_search == 4` before applying threshold of retaining graphs with worst residual absolute values no greater than 3.
- `median_ngraphs_h4_after`: median number of graphs across across 100 replicates when `h_search == 4` retained after applying threshold of retaining graphs with worst residual absolute values no greater than 3.

## popnets

All files were generated from [`summarize_graphs_g2.qmd`](../scripts/summarize_graphs_g2.qmd).

if by paramset, column structure is as follows:

- `net`: name of network
- `h_search`: number of max number of hybridization events searched for (`max_admix` argument to `find_graphs`)
- `trueh`: true number of hybridization events in the network
- `level`: level of true network
- `num_dtrees`: number of displayed trees in true network
- `num_mincyclesize3`: number of cycles with a minimum size of 3 within true network
- `gamma`: gamma value either `large` or `small`. default to `large`, unless a network was simulated under both `large` or `small`, like `g1` / `g4`
- `prop_min_nethw_eq_0`: proportion of replicates that have at least 1 graph with hardwired cluster distance of 0 to true network
- `prop_min_hashmatch`: proportion of replicates that have at least 1 graph that matches the `hash` value of true graph (isomorphic match)
- `prop_min_treehw_eq_0`: proportion of replicates that have at least 1 graph with hardwired cluster distance of 0 to major tree of true graph
- `avgngraphsperrep`: average number of graphs across all replicates
- `prop_max_ndisp_eq_2`: proportion of replicates that have at least 1 graph that displays 2 major trees of the true network (this is not a useful column, was not used to analyze anything and was a holdover from the `g1/g4` version)

if by replicate, column structure  same but with the addition of the `irep` identifier column, plus these replacements after the `gamma` column:

- `ngraphs`: number of graphs retained (varies depending on filtering scheme)
- `min_nethw`: minimum hardwired cluster distance from inferred graphs to true graph
- `atleastonehashmatch`: boolean value whether at least 1 graph that matches the `hash` value of true graph (isomorphic match)
- `min_treehw`: minimum (across replicate) hardwired cluster distance from displayed trees of inferred graph to major tree of true graph
- `max_ndisplayed`: maximum number of displayed trees from inferred graph, with gamma threshold of 0.0, that match near-major trees from true network
(see column `num_dtrees` for correct number of display trees per network)
- `mean_ndisplayed`: mean number of displayed trees from inferred graph, with gamma threshold of 0.0, that match near-major trees from true network
(see see column `num_dtrees` for correct number of display trees per network)
- `relativell`: maximum value of `scorediff` (= score of truenet - score of estimated net), the best score of all inferred graphs

## other files

### Fig. 1 files

both `citations_clean.csv` and `nature_science.csv` have the same column structure:

- `Cites`: number of citations as of Google Scholar query on 11/5/25
- `Authors`: author list
- `Title`: title of article
- `Year`: year of publication (later year taken if merging a preprint and reviewed publication)
- `Source`: journal of publication
- `Publisher`: publisher of publication
- `ArticleURL`: URL of article on Google Scholar
- `CitesURL`: URL of all citations of the article on Google Scholar
- `Type`: type of publication
- `DOI`: DOI, if available
- `CitesPerYear`: `Cites`/`Age`
- `Age`: number of years since publication
- `FullTextURL`: URL of full text if available

### timingdata.csv

This is a file containing the results from [the timingdata.qmd](../scripts/timingdata.qmd) file.
There, I took one parameter combination `fleg_10-indls_118000-genes_100000-biall_0.0-lindist_1.25e-8-subrate` (baseline simulation conditions for g4 graph)
and scraped all the logs containing timing data for graph search:

**timing data column structure**:

- `irep`: replicate number (out of 100 replicates)
- `irun`: run number (ran 50 runs per h_search value in each replicate)
- `h_search`: number of max number of hybridization events searched for (`max_admix` argument to `find_graphs`)
- `time`: amount of minutes before graph search was completed