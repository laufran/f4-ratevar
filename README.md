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
7. `est_batch_g2.jl`: script to run raph search estimation functions from `est_f4_nets.jl` in mass parallelization
      for all other nets besides g1/g4.

## input files

in the `input` folder there are files:

- `fleg_h1_quartets.csv`:
- `fleg_h4_quartets.csv`:
- `fleg-net.jl`:
- `pop-net-properties.csv`:
- `pop-net.jl`:

and directory:

- `pop-net_edgelist`:

## code to analyze & visualize results

in the `scripts` folder there are files:

- `calcwr_g1.jl`:
- `calcwr_g4.jl`:
- `comparegraphs_hwdist0.jl`:
- `correct_ndisp.jl`: 
- `interop_admixtools.jl`:
- `expf2_findgraphs.qmd`
- `expf2_findgraphs.r`
- `expf2_simulated.jl`
- `poissonregression.qmd`
- `probability_012mutations.jl`
- `probability_multiallelic.jl`
- `summarize_f4.qmd`
- `summarize_graphs.qmd`:
- `summarize_graphs_g2.qmd`: 
- `timingdata.qmd`:

## output

the `output` folder contains files:

```
output
├── g4_hwmatches.csv
├── h1_byparamid_filtered.csv
├── h1_byrep_filtered.csv
├── h4_byparamid_filtered.csv
├── h4_byrep_filtered.csv
├── h4_type1_fromwr.csv
└── timingdata.csv
```

### g4_hwmatches.csv

This file contains the graph search "matches" for the `g4` graph,
where they did not match based on `hash` (`truenet == FALSE`),
but do have hardwired cluster distance of 0 (`net_hwdist == 0`).

### h1 & h4 by paramid

### h1 & h4 by rep

### h4_type1_fromwr.csv

### timingdata.csv

This is a file containing the results from [the timingdata.qmd](scripts/timingdata.qmd) file.
There, I took one parameter combination `fleg_10-indls_118000-genes_100000-biall_0.0-lindist_1.25e-8-subrate` (baseline simulation conditions for g4 graph)
and scraped all the logs containing timing data for graph search:

**timing data column structure**:

- `irep`: replicate number (out of 100 replicates)
- `irun`: run number (ran 50 runs per h_search value in each replicate)
- `h_search`: number of max number of hybridization events searched for (`max_admix` argument to `find_graphs`)
- `time`: amount of minutes before graph search was completed
