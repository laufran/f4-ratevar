#=
open julia with `julia --project -p X` w/ X = # of procs
then do : 

```julia
include("scripts/simulate_batch_g2.jl")
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
    include("../input/pop-net.jl")
    using ProgressMeter
end

# nets with different levels, cycles, tree-child vs nongalled, etc.
net_list  = eval.(population_net_names) # defined in pop-net.jl
# obsolete: (g2_l1_c4, g2_l1_c5, g2_l2_c4, g2_l2_c6, g2_l2_ntc, g2_l2_ng4, g2_l2_ng6, g1_sg)
net_names = replace.(string.(population_net_names), "_"=>"-")
# obsolete: ["g2-l1-c4", "g2-l1-c5", "g2-l2-c4", "g2-l2-c6", "g2-l2-ntc", "g2-l2-ng4", "g2-l2-ng6", "g1-sg"]

# 100 replicates per parameter set combination
nreplicates = 100
nrep_threads = 20
nrep_perthread = Int(nreplicates / nrep_threads)

rates_list = collect((n1, 1+nrep_perthread*(irept-1))
    for n1 in 1:length(net_list)
    for irept in 1:nrep_threads)
@show rates_list

p = Progress(length(rates_list), desc="Computing...", barglyphs=BarGlyphs("[=> ]"), showspeed = true)

@sync @distributed for (net, irepstart) in rates_list
    paramset_seed = rand(1:100_000_000)
    @show(paramset_seed)

    simulate_fromnet_tovcf(net_list[net]; nind= 10, lineage_dist= 0.0, 
        sub_rate= 1.25e-8, popsizes= flegpopsizes, nsites= 500, 
        nbiallelic= 100000, ngenes= 118000, ntries= 20, nrep= nrep_perthread, 
        rep_id_start= irepstart, seed= paramset_seed, netname=net_names[net])

    next!(p)
end
