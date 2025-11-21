# pipeline to study accuracy of admixture graph inference from f4 statistics

Code to reproduce the simulation study in: Frankel & Ané (2025) "Low accuracy of complex admixture graph inference from f-statistics".

Preprint available here: <https://www.biorxiv.org/content/10.1101/2025.03.07.642126v2>

## goal and model

Goal: simulate sequence data (SNPs) under a network model reflecting
archaic human population genetics with ILS and lineage rate variation,
to answer these questions:

1. How does f4 perform when there is lineage rate variation?
2. How are the admixture graphs inferred from f4 results affected?
   Specifically, are complex admixture graphs inferred correctly?

## dependencies

- julia packages: see `Project.toml`
- seqgen for the simulation of sequences: Seq-Gen v1.3.4 downloaded from
  [github](git@github.com:rambaut/Seq-Gen.git) at [commit from 2019-08-29](https://github.com/rambaut/Seq-Gen/commit/e8660d73769297d25a88b40d5b828a6fce7784b5).
- R packages:
  - [admixtools](https://github.com/uqrmaie1/admixtools), see [notes/software_install.md](notes/software_install.md#admixtools---done-on-franklin)
- [`vcf2eigenstrat.py`](bin/vcf2eigenstrat.py) Python file, modified from Iain Mathieson's two conversion scripts 
   [here](https://github.com/mathii/gdc/blob/master/eigenstrat2vcf.py) and [here](https://github.com/mathii/gdc/blob/master/gdc.py)

## code to simulate f4-stats + admixture graphs

to simulate data, we used the following scripts in the `scripts` directory:

1. `simulatedata.jl`: functions to simulate from fixed networks (defined in `input/fleg-net.jl` or `input/pop-net.jl`)
      to VCF files. from fixed network, simulates gene trees (with or without lineage rate variation),
      then simulates SNPs from those gene trees, using seq-gen. eventually writes output to VCF.
2. `miscfxns.jl`: contains helper functions used in various scripts throughout simulation and analyzing process.
      examples of functions- converting HybridNetwork objects in Julia to adjacency list format for use with
      admixtools (`net_to_adjlist`), etc.
3. `simulate_batch.jl`: script to run simulation functions from `simulatedata.jl` in mass parallelization for g1/g4.
4. `simulate_batch_g2.jl`: same as above, but for all other nets besides g1/g4. (only ran baseline conditions)
5. `est_f4_nets.jl`: once data (VCF files) are simulated, now we can estimate f4 statistics & estimate/fit admixture
      graphs. contains functions to do so (`calcf4_fromvcf` & `calcgraphs_fromvcf`).
6. `est_batch.jl`: script to run f4/graph search estimation functions from `est_f4_nets.jl` in mass parallelization
      for g1/g4.
7. `est_batch_g2.jl`: script to run graph search estimation functions from `est_f4_nets.jl` in mass parallelization
      for all other nets besides g1/g4.
8. `correct_ndisp.jl`: script to correct a bug in `est_f4_nets.jl`. this file takes our output and recalculates
      the number of displayed true near major trees by an inferred network. the bug has since been corrected,
      so if one were to re-run function `calcgraphs_fromvcf` from `est_f4_nets.jl`, the output would be correct now.
9. `calcwr_g1.jl`: in our original estimation of graphs, we did not calculate the worst residuals of fitting graphs based
      on f4-stats. (it would have been more efficient than doing it afterward, as we have to recalculate f2 blocks. so if
      someone were to do this again, would be better to stick in `calcgraphs_fromvcf` pipeline!) so afterwards, we only took top graphs when inferring 1 or 2 hybridization events (`h_search == {1,2}`) to assess whether users might be
      at risk of overestimating the number of hybridization events when using the worst residual metric.
10. `calcwr_g4.jl`: same as above, but for g4 graph, and `h_search == 4` only.

## input files

in the `input` folder there are files:

- `fleg-net.jl`: defines the g1/g4 graphs (g4 originally coming from [Flegontov et al. 2023](https://elifesciences.org/articles/85492) Fig. 3A)
- `pop-net.jl`: defines all other graphs
- `pop-net-properties.csv`: csv containing properties of all graphs (true number of hyb. events, level, etc.)
- `fleg_h1_quartets.csv`: for every quartet possible in g1 network, contains the number of hybridization events
    before shrinking 2-cycles (h1), after shrinking (h2), and their respective gamma values (g1 for first hybridization
    event, g2 for second, etc...). for more info read `calc_quartet_numhybrids` docstring in [`miscfxns.jl`](scripts/miscfxns.jl)
- `fleg_h4_quartets.csv`: same as above, but for g4

and directory:

- `pop-net_edgelist`: contains a csv for all 18 networks, coding them in edgelist format

## code to analyze & visualize results

in the `scripts` folder there are the files listed below.

main files for data analysis / visualization:

- `summarize_f4.qmd`: markdown doc with analysis of power & type-1 error of f4 statistic across simulation conditions for g1/g4.
- `summarize_graphs.qmd`: markdown doc with analysis of results from `find_graphs` for g1/g4. contains code to concatenate results from the cluster, and generate summaries of data & figures.
- `summarize_graphs_g2.qmd`: same as above, but for all other nets besides g1/g4.

identifiability analyses:

- `expf2_findgraphs.qmd`: markdown document containing calculation of expected f2 values from g4 graph, and fitting f2 values to candidate topologies that could have same score as g4 as part of exploring topological identifiability of g4.
- `expf2_findgraphs.r`: part of analyses in `expf2_findgraphs.qmd`, to fit f2 values to candidate topologies that may have
same score as g4.
- `expf2_simulated.jl`: code to calculate average f2 values, from both observed and simulated data. also contains figure
generation for supplement.
- `interop_admixtools.jl`: contains some conversions of functions in `admixtools`/R to Julia.

supplemental material, etc.

- `probability_012mutations.jl`: script to calculate the marginal probability of 0, 1, 2+ mutations per site,
      across estimated gene trees.
- `probability_multiallelic.jl`: script to estimate the probabilities
      P{site is variable} and P{site is multiallelic | variable} from our simulated sequences.
- `poissonregression.qmd`: markdown document containing calculation of Poisson regression of distance of closest inferred  top graph to true network.
- `comparegraphs_hwdist0.jl`: code to better understand the similarities and differences of the 14 graphs
at hardwired-distance 0 from g4, yet non-isomorphic to g4.
- `timingdata.qmd`: script to take one parameter combination (`fleg_10-indls_118000-genes_100000-biall_0.0-lindist_1.25e-8-subrate`, baseline simulation conditions for g4 graph) and scrape all the logs containing timing data for graph    search (see results below, `timingdata.csv`)

## output

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

Please read the `output` specific readme for more details on column structure, etc. of each file.

### g1-g4

All of these files contain output from [`summarize_graphs.qmd`](scripts/summarize_graphs.qmd).

- `g4_hwmatches.csv`: This file contains the graph search "matches" for the `g4` graph,
where they did not match based on `hash` (`truenet == FALSE`),
but do have hardwired cluster distance of 0 (`net_hwdist == 0`).
- h1 & h4 by rep: .csv file containing results for graph search of g1 (`h1_byrep_filtered.csv`) or g4 (`h4_byrep_filtered.csv`) networks,
as summarized by replicate.
- h1 & h4 by paramid: .csv file containing results for graph search of g1 (`h1_byparamid_filtered.csv`) or g4 (`h4_byparamid_filtered.csv`) networks, as summarized by parameter set.
- `h1_type1_fromwr.csv` and `h4_type1_fromwr.csv`: Summary of results when evaluating type-1 error from using threshold that graphs must have
a worst residual no greater than 3.
Generated from data produced in `calcwr_g1.jl` (results in `h1_type1_fromwr.csv`) and `calcwr_g4.jl`
(results in `h4_type1_fromwr.csv`).

### popnets

All of these files contain results from [`summarize_graphs_g2.qmd`](scripts/summarize_graphs_g2.qmd),
with summary information for the additional 16 population networks beyond g1/g4, but with different filtering schemes.

- `topgraphs_1graph...`: files that summarize graph search results by replicate or parameter set when only retaining the single top scoring graph per replicate
- `topgraphs_3thresh...`: files that summarize graph search results by replicate or parameter set when retaining graphs scoring within 3 of the best score, or 5 min. graphs
- `topgraphs_5graphs...`: files that summarize graph search results by replicate or parameter set when only retaining 5 top scoring graphs per replicate
- `topgraphs_by`: files that summarize graph search results by replicate or parameter set with default filtering scheme (5 minimum graphs, or graphs scoring within 10 of the best score)

### Other files

- `citations_clean.csv`: file containing publication information regarding the query in Fig. 1 (all publications)
- `nature_science.csv`: file containing publication information regarding the query in Fig. 1 (only Nature + Science family publications)
- `timingdata.csv`: This is a file containing the results from [the timingdata.qmd](scripts/timingdata.qmd) file.
There, I took one parameter combination `fleg_10-indls_118000-genes_100000-biall_0.0-lindist_1.25e-8-subrate` (baseline simulation conditions for g4 graph)
and scraped all the logs containing timing data for graph search.