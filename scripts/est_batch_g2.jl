#=
Given simulated data (from `simulate_batch.jl`, resulting in nets through eigenstrat files),
this file will do `admixtools` analyses contained in `calcf4_graphs_fromvcf`:
    - estimate topology of true graph using `qpgraph`, w/ f2_blocks from each replicate as input
    - calculate the best-fitting graphs from f2_blocks for h=0,1, and then the true # of hyb. events
        (minimum of 5 graphs will be saved per h / rep combo, depending on how many graphs
        fall within minimum likelihood score +10)
Output: each replicate will have a folder `graphruns` with files `topgraphsX-XXX.csv` 
    (first X = h_search val, XXX = irun #). should be # of `topgraphsX-XXX.csv` files equivalent to 
    100 runs * # of h_search values per replicate.

to use this script:
open julia with `julia --project -p X` where X = int. # of procs

then do : 
```julia
include("scripts/est_batch_g2.jl")
```
=#

@everywhere begin
    include("est_f4_nets.jl")
    include("../input/pop-net.jl")
end

#get paramset dir names, for use in parallelizing by output directory
output_dirs = sort!(filter!(contains(r"\w+_\d.+"), readdir("output")))

nreplicates = 100
nrep_threads = 50
nrep_perthread = Int(nreplicates / nrep_threads)

rates_list = collect((dir,1+nrep_perthread*(irept-1))
    for dir in output_dirs
    for irept in 1:nrep_threads)
@show rates_list

@sync @distributed for (dir, irepstart) in rates_list
    net  = split(dir, "_")[1]
    net = replace(net, "-" => "_")
    tmpname = (Symbol(net * "_nobot"))
    adjlist = net_to_adjlist(eval(tmpname))
    @show net, irepstart

    if startswith(net, "g3")
        hsearchvals = [0,1,3]
    elseif startswith(net, "g4")
        hsearchvals = [0,1,4]
    elseif startswith(net, "g2") || startswith(net, "g1")
        hsearchvals = [0,1,2]
    else
        @error "unanticipated net name"
    end

    starttime = time()
    #estimate graphs from f2 blocks
    calcgraphs_fromvcf(dir, adjlist, hsearch = hsearchvals, nrep = nrep_perthread, rep_id_start= irepstart)

    endtime = time()
    elapsed = (endtime - starttime) / 3600
    @info "elapsed time: $elapsed hours for repstart: $irepstart, net: $net"
end
