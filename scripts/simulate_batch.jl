#=
open julia with `julia --project -p X` w/ X = # of procs
then do : 

```julia
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
include("scripts/simulate_batch.jl")
```

If running `nbiallelic= 100000, ngenes= 118000` 
 (not `nbiallelic= 10000, ngenes= 11800`), then
 do NOT run all rates at once. On franklin00-02,
 1tb of memory is available and running all rates
 uses ALL OF the memory available and crashes the
 cluster. Request ~75 processors at a time and run
 the rates in batches of 75. (e.g., indexing
 `rates_list on line 48`)
=#

@everywhere begin
    include("simulatedata.jl")
    include("../input/fleg-net.jl")
    using ProgressMeter
end

#fleg or pruned fleg as true net
net_list = (fleg, fleg_pruned)
net_names = ["fleg", "fleg-pruned"]
#1, 2 or 10 individuals per population
nind_list = (1, 2, 10)
#lineage distribution stdev values:
lindist_list = (0.0, 0.15, 0.3, 0.5, 0.7)
#100 replicates per parameter set combination (ie, loci / lindist combo)
nreplicates = 100
nrep_threads = 10
nrep_perthread = Int(nreplicates / nrep_threads)

#I ran one nbi/ngene & mutrate combo on each cluster.
#Uncomment lines below to run diff. combos:
#nbi = 100000; ngene = 118000
nbi = 10000; ngene = 11800
#mutrate = 1.25e-8
mutrate = 1.25e-7

rates_list = collect((n1, i1,s1,1+nrep_perthread*(irept-1))
    for n1 in 1:length(net_list)
    for i1 in nind_list
    for s1 in lindist_list
    for irept in 1:nrep_threads)
@show rates_list

p = Progress(length(rates_list), desc="Computing...", barglyphs=BarGlyphs("[=> ]"), showspeed = true)

#uncomment on next line if running nbi = 100000
    #for each batch, increment slice by 75
    #should be 300 rates in `rates_list` if running one mutrate & nbi/ngene combo at a time
@sync @distributed for (net, nind, lineage_dist, irepstart) in rates_list#[1:75]
    @show(nind)
    @show(lineage_dist)

    paramset_seed = rand(1:100_000_000)
    @show(paramset_seed)

    simulate_fromnet_tovcf(net_list[net]; nind= nind, lineage_dist= lineage_dist, 
        sub_rate= mutrate, popsizes= flegpopsizes, nsites= 500, 
        nbiallelic= nbi, ngenes= ngene, ntries= 20, nrep= nrep_perthread, 
        rep_id_start= irepstart, seed= paramset_seed, netname=net_names[net])

    next!(p)
end
