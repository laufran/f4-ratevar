#=
This file is to calculate worst residual scores for certain datasets:
    -G1 graphs: `fleg-pruned`
    -`h_search` vals: 1, 2

To do so, it first filters the original graph data for the relevant 
parameter sets (`concat_topgraphs.csv` for each replicate). 
This script filters those graphs for the `h_search` value(s) of interest, 
and applies our graph filtering scheme (taking the best scoring graph
by loglikelihood, and the next 4 graphs, or as many graphs fall 
within 10 loglikelihood points). We eventually this apply this filtering
later on (see: `scripts/summarize_graphs.qmd`), but to hasten the process
of calculating worst residuals on only our "best" graphs of interest and 
lower memory needs/help parallelize, we write out these out to 
`topgraphs-acrossiruns.csv` for each replicate.

Then we calculate the worst residuals for these graphs, and write out
a file for each replicate to "output/worstresid/fleg-pruned/$paramset-$irep.csv".
(Again, separate files to help with memory/parallelization.)

After running this on each franklin, the files will need to be concatenated,
zipped, and moved to franklin00 like so:
```{bash}
awk 'FNR==1 && NR!=1{next;}{print}' *.csv > worstresid2.csv
tar -czvf worstresid2.tar.gz worstresid2.csv
scp worstresid2.tar.gz lefrankel@franklin00.stat.wisc.edu:/nobackup2/lefrankel/f4-ratevar/output/
```

to use this script:
open julia with `julia --project -p X` where X = int. # of procs

then do : 
```julia
@everywhere using RCall
@everywhere R"library(admixtools)"
@everywhere R"library(magrittr)"
@everywhere R"library(dplyr)"
include("scripts/calcwr_g1.jl")
```

I couldn't tell you why, but it does not like having the R
packages loaded in an @everywhere block in the file itself.
=#

@everywhere begin
    using DataFrames
    using CSV
    include("miscfxns.jl")

    machine = gethostname()
    if machine == "franklin00"
        paramsets = sort!(filter!(contains(r"fleg-pruned.+118000-genes_100000.+_1.25e-7.+"), readdir("output")))
    else
        paramsets = sort!(filter!(contains(r"fleg-pruned"), readdir("output")))
    end

    reps = range(1,100)
    rates_list = collect((dir,irep)
        for dir in paramsets
        for irep in reps)
end

@sync @distributed for paramset in paramsets
    repdirs = sort!(filter!(contains(r"\d+"), readdir("output/$paramset")))

    for rep in repdirs
        @show "on paramset: $paramset rep: $rep"
        
        #load in CSV
        irep_df = DataFrame(CSV.File("output/$paramset/$rep/concat_topgraphs.csv"))

        #filter out h_search == 0
        filter!(:h_search => !isequal(0), irep_df)        
        #filter out true graph
        filter!(:truenet => isequal(false), irep_df)
        #group by h_search
        groups_byhsearch = groupby(irep_df, [:h_search])
        #convert subdataframes to dataframes so I can use modifying functions below
        groups_byhsearch = [DataFrame(groups_byhsearch[1]), DataFrame(groups_byhsearch[2])]
        #initialize df
        irep_df = DataFrame()

        for group in groups_byhsearch
            #first sort on score
            sort!(group, :score)
            #keep only unique hashes
            unique!(group, :hash)
            #get min score
            min_score = group.score[1] #first row will be smallest score
            threshold_score = maximum([min_score + 10, group.score[5]])
            #filter within that threshold score
            group = group[group.score .<= threshold_score, :]
            #push to df
            append!(irep_df, group)
        end

        #rewrite filtered csv out
        CSV.write("output/$paramset/$rep/topgraphs-acrossiruns.csv", irep_df)
    end
end

@sync @distributed for (paramset, irep) in rates_list
    starttime = time()

    #check whether files exist first as error handling
    if isfile("output/worstresid/$paramset-$irep.csv"); continue; end

    repstr = string(irep, pad = 3)
    #calculate f2 blocks from alignment
    rep_dir = "output/" * paramset * "/" * repstr
    prefix = rep_dir * "/seq"
    R"f2_blocks = f2_from_geno($prefix, blgsize = 10000, adjust_pseudohaploid = TRUE, verbose = FALSE)"
    
    #initialize wr list per irep
    wr_list = []

    #load in CSV
    irep_df = DataFrame(CSV.File("output/$paramset/$repstr/topgraphs-acrossiruns.csv"))

    #convert newick cols to igraph
    newicks = irep_df.newick
    for tree in newicks
        #convert newick to adj list
        edges = net_to_adjlist(tree)
        @rput edges
        #calculate worst residual using code from Rob Maier github issue
        R"wr = qpgraph(f2_blocks, edges, return_fstats = T)$f4 %>%
            slice_max(abs(z), with_ties = F) %>% pull(z)"
        @rget wr
        #push wr to list
        push!(wr_list, wr)
    end

    #insert list as col in df
    insertcols!(irep_df, 3, :worstresid => wr_list)

    #update writing df as CSV
    CSV.write("output/worstresid/fleg-pruned/$paramset-$irep.csv", irep_df)

    endtime = time()
    elapsed = (endtime - starttime) / 3600
    @show "paramset: $paramset rep: $irep, elapsed time: $elapsed hours"
end