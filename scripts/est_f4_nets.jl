#=
Code for analyzing data using admixtools:
1. `calcf4_fromvcf`: outputs `f4_d.csv` per replicate
2. `calcgraphs_fromvcf`: outputs separate graph search files per run & `search_h` value in 
    `graphreps` directory, within each replicate directory.

Helper functions:
1. `ndisplayed_nearmajortrees`

This function is called & parallelized by `est_batch.jl`.
=#

include("../input/fleg-net.jl")
include("miscfxns.jl")
using PhyloNetworks
using RCall
using StatsBase
using Dates
using CSV
using DataFrames

#load R deps
#turning warnings off about masking diff. fxns for 
  #`magrittr` and `dplyr` because it's annoying when
  #running `est_batch.jl` on 100 procs
R"library(admixtools)"
R"library(stringr)"
R"suppressWarnings(library(dplyr))"
R"library(readr)"
R"library(rlang)"
R"suppressWarnings(library(magrittr))"

"""
    ndisplayed_nearmajortrees(net::HybridNetwork, trees2find=nearmajortree, nmt=num_nearmajortrees)

Number of trees to find in the list `trees2find` that are displayed by `net` (in HybridNetwork format).
Arguments and default values (global variables, to avoid being re-calculated):
- `trees2find = nearmajortree`
- `nmt = num_nearmajortrees`: length of `trees2find`
"""
function ndisplayed_nearmajortrees(net, trees2find=nearmajortree, nmt=num_nearmajortrees)
  dtrees = removedegree2nodes!.(displayedTrees(net, 0.0))
  notfound = collect(eachindex(trees2find))
  for tree in dtrees
    isempty(notfound) && break
    for (ii,i) in enumerate(notfound)
      if hardwiredClusterDistance(tree, trees2find[i], false) == 0
        deleteat!(notfound, ii) # no need to look for trees2find[i] in later displayed trees
        break # 'tree' cannot match other near-major trees
      end
    end
  end
  return nmt - length(notfound)
end

"""
    function calcf4_fromvcf(path::AbstractString;
      nrep::Int, rep_id_start::Int)

Given simulated eigenstrat data (from `simulate_fromnet_tovcf`, should be 3 files
called seq.snp, seq.ind, and seq.geno) for each replicate in a output directory located
at `path`, this function will do `admixtools` analyses:
  - calculate f4 / D statistics for every four-taxon combination. `adjust_pseudohaploid` is set
    to TRUE, as we simulated haploid individuals and only have 0s or 2s in the eigenstrat files
    (no 1s) --> f4 then treats each individual as *haploid*, not diploid.
Output: each replicate will have 1 csv:
  - `f4_d.csv`: contains cols `pop1,pop2,pop3,pop4,f4,f4_se,f4_Z,f4_P,d,d_se,d_Z,d_P,irep,lindist,nind`.
      f4 & d statistics, their standard error, Z- and P-vals are calculated for every quartet. 
      paramset info (`nind` and `lindist`) as well as rep # (`irep`) are appended.

For parallelization purposes, specify number of replicates to analyze `nrep` and
replicate to start on `rep_id_start`. For example, if no parallelization (not recommended)
and there are 100 replicates produced by `simulate_fromnet_tovcf`, specify `nrep` = 100
and `rep_id_start` = 1.
"""

function calcf4_fromvcf(path::AbstractString;
            nrep::Int, rep_id_start::Int)

  dirstring = split(path, "_")
  nind =      parse(Int, split(dirstring[2], "-")[1])
  lindist =   parse(Float64, split(dirstring[5], "-")[1])
  
  irep = rep_id_start
  rep_id_end = rep_id_start + nrep
  while irep < rep_id_end
    repstring = string(irep, pad = 3)

    if isfile("output/$path/$repstring/f4_d.csv")
      irep += 1
      continue
    end

    f4_d_df = DataFrame(pop1= String[], pop2= String[], pop3= String[], pop4= String[], 
                    f4= Float64[], f4_se= Float64[], f4_Z= Float64[], f4_P= Float64[],
                    d=  Float64[], d_se=  Float64[], d_Z=  Float64[], d_P=  Float64[],
                    irep= Int64[], lindist= Float64[], nind= Int64[])

    #get f2 blocks
    rep_dir = "output/" * path * "/" * repstring
    prefix = rep_dir * "/seq"
    #blgsize = 10000 regardless of number of SNPs, because we artificially
    #incremented each SNP in the VCF by 10000bp
    R"f2_blocks = f2_from_geno($prefix, blgsize = 10000, adjust_pseudohaploid = TRUE)"
    R"if(length(f2_blocks) < 9990) warning(str_glue('Have less SNPs than intended (only ', length(f2_blocks),' SNPs), check input eigenstrat files.'))"
    
    #calculate f4
    R"f4_alltaxa = f4(f2_blocks, f4mode = TRUE)"
    R"f4_alltaxa = f4_alltaxa %>% rename(f4 = est, f4_se = se, f4_Z = z, f4_P = p)"

    #calculate D
    R"d_alltaxa = f4(f2_blocks, f4mode = FALSE)"
    R"d_alltaxa = d_alltaxa %>% rename(d = est, d_se = se, d_Z = z, d_P = p)"

    #merge dfs- only if the column identies of d & f4 match
    R"stopifnot(d_alltaxa[,1:4]==f4_alltaxa[,1:4])"
    R"f4_d_alltaxa = cbind(f4_alltaxa, d_alltaxa[!names(d_alltaxa) %in% names(f4_alltaxa)])"
    R"f4_d_alltaxa$irep = $irep" #insert column with rep
    R"f4_d_alltaxa$lindist = $lindist" #insert column with lindist
    R"f4_d_alltaxa$nind = $nind" #insert column with nind
    #get f4_d_alltaxa R --> Julia, append to f4_d_df Julia obj
    @rget f4_d_alltaxa
    f4_d_df= append!(f4_d_df, f4_d_alltaxa)

    CSV.write("$rep_dir/f4_d.csv", f4_d_df)

    irep += 1
  end #for each rep
end

"""
    calcgraphs_fromvcf(path::AbstractString, net_adjlist::DataFrame, nruns::Int=50;
      hsearch::Vector{Int}, nrep::Int, rep_id_start::Int)

Given simulated eigenstrat data (from `simulate_fromnet_tovcf`, should be 3 files
called seq.snp, seq.ind, and seq.geno) for each replicate in a output directory located
at `path`, an adjacency list for the true network `net_adjlist`, this function
will do `admixtools` analyses:
- estimate edge parameters (lengths, gammas) on the true graph topology using
  `qpgraph`, with f2_blocks from each replicate as input
- search for best-fitting graph topologies (and edge parameters) from f2_blocks
  for hs in `hsearch` vector `nruns` times (runs), 50 by default
  A minimum of 5 graphs will be saved per h / rep combo. More can be saved:
  all graphs whose score falls within minimum likelihood score +10.
  default: 5*50 runs -> 250 minimum graphs per replicate

Output: each replicate will have 1 log:
  - graphsearch_log.txt`: contains the graph search arguments, seed per `find_graphs` run 
      & elapsed time (should be 50 (or `nruns`) runs per replicate & h_search val)
  and one .csv per run:
  - `topgraphsX-XXX.csv` (first X = h val, XXX = run #): contains cols:
    - `newick`: Newick string for each tree/net
    - `score`: log-likelihood for given tree/net
    - `hash`: given from `admixtools`, to (maybe) be used to see if topologies are isomorphic
    - `truenet`: boolean (if TRUE: true network est. from `qpgraph`, if FALSE: estimated from 
      `find_graphs`)
    - `net_hwdist`: hardwired cluster distance from true net (0 if net is true net)
    - `tree_hwdist`: hardwired cluster distance from major tree of true graph (0 if net is true net)
    - `scorediff`: `truenet`- "estimated score"'s score (0 if net is true net) 
        this is the loglikelihood of the estimated graph minus the loglikelihood
        of the true topology. higher is better!
    - `h_search`: number of hybridization events searched for (`find_graphs` arguments
        `numadmix` / `max_admix`). if true net, fleg will have h_search of 4, fleg-pruned 1
    - `h_est`: number of hybridization events estimated. if true net, fleg will have h_search
        of 4, fleg-pruned 1
    - `ndisplayed`: number of displayed trees from inferred graph, with gamma threshold of 0.0,
        that match either of 2 near-major trees from true network (gamma > 0.45). values should
        not exceed 2. (0 = bad, 2 = best)
    - `irep`: replicate # graph simulated from
    - `irun`: run # graph simulated from (check `graphsearch_log.txt` to see seed per run)
    - `lindist`: standard deviation for lognormal distribution with mean one
    - `nind`: number of individuals simulated per population

For parallelization purposes, specify number of replicates to analyze `nrep` and
replicate to start on `rep_id_start`. For example, if no parallelization (not recommended)
and there are 100 replicates produced by `simulate_fromnet_tovcf`, specify `nrep` = 100
and `rep_id_start` = 1.
"""

function calcgraphs_fromvcf(path::AbstractString, net_adjlist::DataFrame, nruns::Int=50;
  hsearch::Vector{Int}, nrep::Int, rep_id_start::Int)

  length(unique(hsearch)) == length(hsearch) || error("need all unique values in `hsearch`")
  dirstring = split(path, "_")
  nind =      parse(Int, split(dirstring[2], "-")[1])
  lindist =   parse(Float64, split(dirstring[5], "-")[1])
  truenet =   split(path, "_")[1]
  if truenet == "fleg"
    truenet = fleg_nobott
  elseif truenet == "fleg-pruned"
    truenet = fleg_pruned_nobott
  else
      @error "need to specify adjacency list"
  end
  
  nearmajortree = displayedTrees(truenet, 0.45)
  num_nearmajortrees = length(nearmajortree)

  #initialize number of runs for each graph search
  numruns = collect(1:nruns)

  irep = rep_id_start
  rep_id_end = rep_id_start + nrep
  while irep < rep_id_end
    repstring = string(irep, pad = 3)

    if isdir("output/$path/$repstring/graphreps") 
      resultcsvs = sort!(filter!(contains(r"topgraphs"), readdir("output/$path/$repstring/graphreps")))
      #length of resultcsvs may be greater than anticipated, if ran different h_search vals in the past with
      #if length(resultcsvs) >= (nruns * length(hsearch))
      #  irep += 1
      #  continue #skip onto next irep, all iruns are completed & written to csv
      #end
    else #is no `graphreps` directory, need to make one and get "true" graph and estimate all other reps / iruns
      mkdir("output/$path/$repstring/graphreps")
    end

    #get f2 blocks
    rep_dir = "output/" * path * "/" * repstring
    prefix = rep_dir * "/seq"
    #blgsize = 10000 regardless of number of SNPs, because we artificially
    #incremented each SNP in the VCF by 10000bp
    R"f2_blocks = f2_from_geno($prefix, blgsize = 10000, adjust_pseudohaploid = TRUE)"
    R"if(length(f2_blocks) < 9990) warning(str_glue('Have less SNPs than intended (only ', length(f2_blocks),' SNPs), check input eigenstrat files.'))"

    #do qpgraph on true graph adj list using f2_blocks
    R"net_qpgraph = qpgraph(f2_blocks, $net_adjlist)"
    #get hash of true graph from igraph
    R"hash = graph_hash(edges_to_igraph(net_qpgraph$edges))"
    @rget net_qpgraph hash
    net_edges = net_qpgraph[:edges]
    net_score = net_qpgraph[:score]
    #= true graph: always same topology (the one used to simulate data)
                   because it's the topology that was given to qpGraph.
    but: edge lengths and gamma's best fitted the simulated data,
         varies across replicates =#
    truegraph = edgelist_to_net(net_edges[!, "from"], net_edges[!, "to"], 
                  net_edges[!, "type"], net_edges[!, "weight"])
    net_h = truegraph.numHybrids
    truetree = nothing
    try truetree = majorTree(truegraph) # does not include the degree-2 bottleneck node bc from qpGraph
    catch e
      @error """cannot extract major tree from true graph with estimated params, from:
      $rep_dir
      problematic truegraph:
      $(writeTopology(truegraph))
      """
    end
    #= this tree depends on the replicate through the estimated gamma's:
       keeps hybrid edges with gamma > 0.5. It may NOT be the true major tree
       from the true network. But okay, especially for the flegontov network with
       partner hybrid edges with gamma very close to 0.5. =#
    truegraph_string = writeTopology(truegraph, internallabel=false)

    #start logging (per replicate, with info from each run in there)
    #having one log per replicate (not per rep/h_search) might be a problem
      #if i want to run the same replicate but different hs at the same time
    logio = open("$rep_dir/graphsearch_log.txt", "a+") #append if file already exists
    datetime = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    write(logio, "$datetime\n")
    write(logio, """
    stop_gen = 10000
    stop_gen2 = 60
    plusminus_generations = 100
    numstart = 100
    numgraphs = 25
    ------
    """)

    #run find_graphs on f2 blocks, saving min. top 5 scoring graphs per run
    # (more than 5 if graphs fall in min+10 ll score threshold)
    for h in hsearch
      if @isdefined resultcsvs
        stringh = string(h)
        h_csvs = filter(contains(Regex("topgraphs$(stringh)-")), resultcsvs)
        if isempty(h_csvs) #no files exist for this h
          numruns = collect(1:nruns)
        else
          #extract `irun` numbers from existing files
          irun_numbers = Int[]
          for csvfile in h_csvs
            strsplit = match(r"(\d+)-(\d+)\.csv", csvfile)
            push!(irun_numbers, parse(Int, strsplit[2]))
          end
    
          #sort and detect missing `irun` numbers
          irun_numbers = sort(irun_numbers)
          expected_runs = collect(1:nruns)
          numruns = setdiff(expected_runs, irun_numbers) #find missing iruns
    
          if isempty(numruns)
            continue #all runs for this h are complete, go to the next `h`
          end
        end
      else #no files exist previously (resultcsvs isn't defined)
        numruns = collect(1:nruns)
      end

      #generate seeds for all runs per h
      seeds = sample(1:100_000, nruns; replace=false)

      for irun in numruns
        graphs_df = DataFrame(newick= String[], score= Float64[], hash= String[],
          truenet= Bool[], net_hwdist= Float64[], tree_hwdist= Float64[], 
          scorediff= Float64[], h_search= Int[], h_est= Int[], ndisplayed = Int[],
          irep=Int64[], irun=Int64[], lindist=Float64[], nind= Int[])
        allowmissing!(graphs_df, :ndisplayed)
        push!(graphs_df, (truegraph_string, net_qpgraph[:score], hash,
          true, 0.0, 0.0, 0.0, net_h, net_h, missing, 
          #truenet, net_hwdist, tree_hwdist, scorediff, h_search, h_est, ndisplayed
          irep, 0, lindist, nind)) #0 for `irun`

        starttime = time()
        #index seed
        irunseed = seeds[irun]
        #search for graphs
        #if any of these search params change, change hardcoded log lines above
        R"graphs = find_graphs(f2_blocks, outpop = 'chimp', stop_gen = 10000, stop_gen2 = 60, plusminus_generations = 100, numstart = 100, numgraphs = 25, numadmix = $h, seed = $irunseed, return_fstats = TRUE) %>% arrange(score)"
        #get threshold: either min+10 or fifth smallest score, whichever is larger
        R"l_threshold = max(graphs$score[[1]]+10, graphs$score[[5]])"
        #filter for graphs with scores in that threshold
        R"top_graphs = graphs %>% filter(score <= l_threshold)"
        R"top_graphs = subset(top_graphs, select = -c(graph, generation, mutation, lasthash))"
        R"top_graphs$h_search = $h" #insert column with h_search
        R"top_graphs$irep = $irep" #insert column with irep
        R"top_graphs$irun = $irun" #insert column with irun
        R"top_graphs$lindist = $lindist" #insert column with lindist
        R"top_graphs$nind = $nind" #insert column with nind
        R"top_graphs = distinct(top_graphs)" #remove duplicate rows from each irun

        @rget top_graphs 
        top_graphs_edges = top_graphs[!, :"edges"]
        top_graphs_newick = []
        top_graphs_nethwdist = []
        top_graphs_treehwdist = [] 
        top_graphs_hest = []
        top_graphs_ndisplayed = []
        for edgelist in top_graphs_edges
          net = edgelist_to_net(edgelist[!, "from"], edgelist[!, "to"], 
                    edgelist[!, "type"], edgelist[!, "weight"])
          #= 
          originally we did not do `removedegree2nodes!(net, true)`, which led to a bug
          in calculating the hardwired cluster distance from the true net to estimated.
          `find_graphs`` sometimes adds degree-2 nodes because of its restriction to fix hybrid
          edge lengths to 0. These nodes should be suppressed before calculating the distance 
          to the true network, with removedegree2nodes!(net, true).

          i've corrected this below in case this needs to be re-run eventually.
          but the `correct_nethwdist.jl` script is a patch on our buggy version.
          =#
          removedegree2nodes!(net, true) # keeproot=true
          push!(top_graphs_newick, writeTopology(net, internallabel=false))
          push!(top_graphs_hest, net.numHybrids)
          push!(top_graphs_ndisplayed, ndisplayed_nearmajortrees(net, nearmajortree, num_nearmajortrees))
          push!(top_graphs_nethwdist, hardwiredClusterDistance(truegraph, net, true))

          #=
          originally we had `displayedTrees(net, 0.1)` and `nearmajortree[1]` as bugs.
          instead we want to display *all* trees from the estimated network
          with no gamma threshold (= setting to 0.0).
          and erroneously we were looking at the "first" major tree, but to match the viz
          we have in figures of turning one minor gamma to display as major, we should look
          at the second major tree `nearmajortree[2]`.

          i've corrected this below in case this needs to be re-run eventually.
          but the `correct_treehwdist.jl` script is a patch on our buggy version.
          =#

          minhwd = Inf
          if isnothing(truetree)
            @warn "cannot calculate min displayed tree distance to missing truetree"
          else
            treesdisplayedingraph = HybridNetwork[]
            treesdisplayedingraph = displayedTrees(net, 0.0)
            for tree in treesdisplayedingraph
              #checking against just the first near major tree
              hwtreedist = hardwiredClusterDistance(nearmajortree[2], tree, true) # true: rooted at same outgroup
              if hwtreedist < minhwd
                minhwd = hwtreedist
              end
              minhwd == 0 && break # do not look at other displayed trees
            end
          end
          push!(top_graphs_treehwdist, minhwd)
        end

        insertcols!(top_graphs, 1, :newick => top_graphs_newick)
        select!(top_graphs, Not(:edges))
        insertcols!(top_graphs, 5, :truenet => false)
        insertcols!(top_graphs, 6, :net_hwdist => top_graphs_nethwdist)
        insertcols!(top_graphs, 7, :tree_hwdist => top_graphs_treehwdist)
        scorediff_list = net_score .- top_graphs[!, :score]
        insertcols!(top_graphs, 8, :scorediff => scorediff_list)
        #skip a number so this comes after h_search
        insertcols!(top_graphs, 10, :h_est => top_graphs_hest)
        insertcols!(top_graphs, 11, :ndisplayed => top_graphs_ndisplayed)
          
        #add top_graphs to irun specific df
        append!(graphs_df, top_graphs)

        irunstring = string(irun, pad = 3)
        CSV.write("$rep_dir/graphreps/topgraphs$h-$irunstring.csv", graphs_df)

        endtime = time()
        elapsed = (endtime - starttime) / 60 #in minutes

        write(logio, "irun: $irunstring, h_search: $h, seed: $irunseed, time: $elapsed minutes\n")
        flush(logio)
      end #for each graph search run
    end #for each graph search h

    close(logio)
    irep += 1
  end #for each rep
end