#=
Given simulated data (from `simulate_batch.jl`, resulting in nets through eigenstrat files),
this file will do `admixtools` analyses contained in `calcf4_graphs_fromvcf`:
    - calculate f4 / D statistics for every four-taxon combination
    - estimate topology of true graph using `qpgraph`, w/ f2_blocks from each replicate as input
    - calculate the best-fitting graphs from f2_blocks from h=0 -> 4 (true # of hyb. events)
        (minimum of 5 graphs will be saved per h / rep combo, depending on how many graphs
        fall within minimum likelihood score +10)
Output: each replicate will have csv `f4_d.csv` and a folder `graphruns` with files `topgraphsX-XXX.csv` 
    (first X = h_search val, XXX = irun #). should be # of `topgraphsX-XXX.csv` files equivalent to 
    100 runs * # of h_search values per replicate.

to use this script:
open julia with `julia --project -p X` where X = int. # of procs

then do : 
```julia
include("scripts/est_batch.jl")
```
=#

@everywhere begin
    include("est_f4_nets.jl")
    include("../input/fleg-net.jl")
end

#get paramset dir names, for use in parallelizing by output directory
output_dirs = sort!(filter!(contains(r"\w+_\d.+"), readdir("output")))
#special regex for franklin00
#output_dirs = sort!(filter!(contains(r".+118000-genes_100000.+_1.25e-7.+"), readdir("output")))

nreplicates = 100
nrep_threads = 10
nrep_perthread = Int(nreplicates / nrep_threads)

rates_list = collect((dir,1+nrep_perthread*(irept-1))
    for dir in output_dirs
    for irept in 1:nrep_threads)
@show rates_list

@sync @distributed for (dir, irepstart) in rates_list
    dirstring = split(dir, "_")
    nind =      parse(Int, split(dirstring[2], "-")[1])
    lindist =   parse(Float64, split(dirstring[5], "-")[1])
    @show nind, lindist, irepstart

    starttime = time()

    #calculate f4 & d results into 1 csv
    calcf4_fromvcf(dir, nrep = nrep_perthread, rep_id_start= irepstart)

    #specify adjacency lists & h_search values per graph
    if startswith(dir, "fleg_")
        adjlist = net_to_adjlist(fleg_nobott)
        hvals = [0,1,4]
    elseif startswith(dir, "fleg-pruned_")
        adjlist = net_to_adjlist(fleg_pruned_nobott)
        hvals = [0,1,2]
    else
        @error "need to specify adjacency list"
    end

    #estimate graphs from f2 blocks
    calcgraphs_fromvcf(dir, adjlist, hsearch = hvals, nrep = nrep_perthread, rep_id_start= irepstart)

    endtime = time()
    elapsed = (endtime - starttime) / 3600
    @info "elapsed time: $elapsed hours for repstart: $irepstart, nind: $nind & lindist: $lindist"
end