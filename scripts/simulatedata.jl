#= 
code to simulate data:
1. `simulatenetwork`: simulate networks using SiPhyNetwork (via RCall)
2. `scale_network`: scale network to control the amount of ILS
3. `simulategenetrees`: simulate gene trees within a network, using PhyloCoalSimulations, 
    and simulate lineage-rate variation if applicable. if given a simulated network,
    scale # of generations from the ingroup crown to the ingroup taxa to be the same 
    across all generated networks
4. `sim_snps`: simulates SNPs from gene trees, using seqgen

code to analyze simulated data & helper functions:
    - `save_lineagerates`: write lineage rates simulated in `simulategenetrees`
    - `lognormal_meanone`: generate a lognormal distribution with mean of 1
    - `outgrouplabels`: determine outgroups vs. ingroups 
    - `clustersize`: get size of the edge's hardwired cluster (required for `outgrouplabels`)
    - `make_indfile`: given list of individuals, makes a map file of populations for use in converting vcf -> eigenstrat
    - `parameterset_directory`: generate directory string for output files given a set of parameters 
    - `readmultigts_seqgen`: reads gene trees when written to seq-gen format
    - `simulate_fromnet_tovcf`: run through simulation pipeline from a starting network to eigenstrat files

import the file to use its functions:
`include("simulatedata.jl")`
=#

using Random
using PhyloNetworks
using PhyloCoalSimulations
using RCall
using Distributions
using CSV
using DataFrames
using Dates
using Logging
using BioSequences
using FASTX
include("miscfxns.jl")

# below: add with full path to executable if "seq-gen" is not on the path
host = chomp(read(`hostname`, String))
if occursin("Lauren", host)
    seqgen = "/Users/laurenfrankel/Desktop/f4-ratevar/bin/Seq-Gen-1.3.4/source/seq-gen"
elseif occursin("franklin", host)
    seqgen = "/u/l/e/lefrankel/bin/Seq-Gen-1.3.4/source/seq-gen"
elseif occursin("darwin", host)
    seqgen = "/u/l/e/lefrankel/bin/Seq-Gen-1.3.4/source/seq-gen"
elseif occursin("cecile", host)
    seqgen = chomp(read(`which seq-gen`, String))
else
    @warn "Need to specify Seq-Gen location"
end
ispath(seqgen) || error("can't find seqgen at $seqgen")

# R"install.packages('devtools')"
# R"require(devtools)";
# R"devtools::install_github('jjustison/SiPhyNetwork')";
# R"devtools::install_github('uqrmaie1/admixtools')";
# use 10 as the precision for the disribution of γ values
# setting genetic distance function using example values
R"require(ape, quietly=T)";
R"require(SiPhyNetwork, quietly=T)";
R"require(admixtools, quietly=T)";
R"inheritance_fxn <- make.beta.draw(2, 2)";
R"gen_dist_fxn = make.stepwise(probs = c(1,0),distances = c(0.75,Inf))";
hybridproportions = [0.5,0.25,0.25]  # generative, degenerative, neutral
@rput hybridproportions

using Pkg
Pkg.status(["PhyloNetworks","RCall","PhyloCoalSimulations"]);
R"rversion <- R.Version()$version.string";
@rget(rversion);
R"spnversion <- sessionInfo(package='SiPhyNetwork')$otherPkgs$SiPhyNetwork$Version";
R"spnsha     <- sessionInfo(package='SiPhyNetwork')$otherPkgs$SiPhyNetwork$RemoteSha";
@rget(spnversion);
@rget(spnsha);
R"admixversion <- sessionInfo(package='admixtools')$otherPkgs$admixtools$Version";
R"admixsha     <- sessionInfo(package='admixtools')$otherPkgs$admixtools$RemoteSha";
@rget(admixversion);
@rget(admixsha);
python_version = ""
try
    global python_version = readchomp(`python --version`)
catch _
    @warn "'python' not in your path"
end
print("""julia version: $VERSION
      $rversion
      SiPhyNetwork version: v$spnversion
      SiPhyNetwork commit sha: $spnsha
      admixtools version: v$admixversion
      admixtools commit sha: $admixsha
      Python version: $python_version
      """)

"""
    simulatenetwork(ntaxa::Int, nsims::Int)

calls `simulatenetwork(ntaxa)` many times, to simulate `nsims` networks
and return them as a tuple with vector of `HybridNetworks` and a vector of 
levels of the networks. The keyword arguments are passed, except for the seed. 
"""
function simulatenetwork(ntaxa::Int, nsims::Int; seed=nothing, kwargs...)
    nets = HybridNetwork[]
    levels = Int[]
    if isnothing(seed)
        seed = rand(1:10000)
        @info "seed for the generation of random seeds: $seed"
    end
    Random.seed!(seed)
    s = rand(1:100_000, nsims) # each s value will be used as a seed within R for 1 network
    for i in 1:nsims
        n, l = simulatenetwork(ntaxa; seed = s[i], kwargs...)
        push!(nets, n)
        push!(levels, l)
    end
    return nets, levels
end

"""
    simulatenetwork(ntaxa::Int, file::AbstractString; seed=nothing,
        lambda=1, mu=0.2, nu=0.2,
        shrink=true, hmin=1, hmax=ntaxa/2, levelmax=hmax)

Simulate 1 network using SiPhyNetwork repeatedly until the network
satisfies the filters, append it to `file` and return it as a tuple of 
`HybridNetwork` and level of network.

The filters are:
- exclude networks with fewer than `hmin` reticulations (trees by default)
- exclude networks with more than `hmax` reticulations
- exclude networks of level higher than `levelmax`

Before checking the level, there is the option to `shrink` 2-cycles and 3-cycles
(which are not detectable by the ABBA-BABA test or by SNaQ).

By default, the seed is not controlled. Otherwise, it is passed to R to seed
the random number generator for the SiPhyNetwork simulation.

external variables (which could be turned into options):
- proportion of the 3 kinds of reticulations: `hybridproportions`
- distribution of γs: R function `inheritance_fxn`
- probability that a proposed hybridization is successful, given the
  average genetic distance between the 2 taxa: R function `gen_dist_fxn`

# example

```repl
julia> network, level = simulatenetwork(7, file="test.phy", seed=121, verbose=false)

julia> network
HybridNetwork, Rooted Network
15 edges
15 nodes: 7 tips, 1 hybrid nodes, 7 internal tree nodes.
tip labels: t7, t10, t9, t12, ...
((((t7:0.018,t10:0.018)I1:0.358,#H15:0.0::0.155)I2:0.195,(t9:0.376,(t12:0.376)#H15:0.0::0.845)I3:0.195)I4:0.984,((t5:0.011,t6:0.011)I5:0.374,t1:0.385)I6:1.17)I7;

julia> level
1
```

R example with dangling nodes, errors with write.net etc:
```r
require(SiPhyNetwork)
inheritance_fxn = make.beta.draw(10,10)
hybridproportions = c(0.5,0.25,0.25)
set.seed(5282)
net <- sim.bdh.taxa.ssa(n=17, numbsim=1,
    lambda=1,mu=0.2, nu=0.8, hybprops=hybridproportions,
    hyb.inher.fxn = inheritance_fxn)[[1]]
write.net(net, digits=1)

library(SiPhyNetwork)
gen_dist_fxn = make.stepwise(probs = c(1,0),distances = c(0.75,Inf))
set.seed(16)
net <- sim.bdh.taxa.ssa(n=4, numbsim=1,
    lambda=1, mu=0.2, nu=0.3, hybprops=c(0.5,0.25,0.25),
    hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = make.beta.draw(2,2), complete=FALSE)[[1]]
write.net(net, tol=1e-12)

# how to find a seed for a minimal working example:
for (i in 1:1000){
  cat("i =",i,"\n")
  set.seed(i)
  net <- sim.bdh.taxa.ssa(n=4, numbsim=1,
    lambda=1, mu=0.2, nu=0.3, hybprops=c(0.5,0.25,0.25),
    hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = make.beta.draw(2,2), complete=FALSE)[[1]]
  if (length(net) < 2){
    cat("went extinct\n")
    next
  }
  write.net(net, tol=1e-12)
}
```
"""
function simulatenetwork(ntaxa::Int; file=nothing, seed=nothing,
        lambda=1, mu=0.2, nu=0.3,
        shrink=true, hmin=1, hmax=ntaxa/2, levelmax=hmax, verbose=true)
    levelmax>0 || error("a non-tree network must have level 1 or more")
    hmin <= hmax || error("hmin may not be higher than hmax")
    isnothing(seed) || R"set.seed"(seed)
    net = nothing
    level = nothing
    for i in 1:1000
        i += 1
        R"""
        net <- sim.bdh.taxa.ssa(n=$ntaxa, numbsim=1,
            lambda=$lambda, mu=$mu, nu=$nu, hybprops=hybridproportions,
            hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = inheritance_fxn, complete = FALSE)[[1]]
        """
        # continue to next iteration if the lineage went completely extinct
        rcopy(R"length(net) > 1") || continue
        # convert R network to Julia network
        net_string = rcopy(R"write.net(net, tol=1e-12)")
        net = readTopology(net_string)
        shrink && PhyloNetworks.shrink3cycles!(net)
        h = net.numHybrids
        (h>=hmin && h<=hmax) || continue # to next iteration if h outside range
        bi = PhyloNetworks.blobInfo(net)[3]
        level = maximum(length.(bi))
        level <= levelmax || continue
        verbose && @info "$i iterations, level=$level"
        break
    end
    PhyloNetworks.nameinternalnodes!(net, "I") # to get be able to map internal node names in gene trees later
    isnothing(file) || writeTopology(net, file; append=true)
    return net, level
end

"""
    scale_network(net::HybridNetwork, file, med_coal_unit::AbstractFloat)

Given a network `net`, median coalescent unit to scale all branches to `med_coal_unit`,
and file path to write network to `file`, `scale_network` scales all branch lengths
in `net` by float `scaler` and writes to `file`.
The goal is to control the amount of ILS.
The amount of substitutions is controlled later.
"""
function scale_network(net::HybridNetwork, file, med_coal_unit::AbstractFloat)
    open(file, "w") do fout
        scaler = med_coal_unit / median([e.length for e in net.edge])
        for e in net.edge
            e.length *= scaler
        end
        writeTopology(net, fout)
    end # close file cleanly
    return net
end

"""
    simulategenetrees_coalescentunits(net, nloci, genefile;
        lineage_distribution, nodemapping=false, nsites=1)

**Not used so far. May be used when networks are simulated at random.
Delete if we end up not using this at all.**

Simulate `nloci` gene trees (1 individual/species) along `net`
using PhyloCoalSimulations, assuming that edge lengths are in coalescent units.
The gene trees with branch lengths in coalescent units are written to file
"`genefile`_coal_unit". These gene trees have extra degree-2 nodes to map
them into the species network.

The model for rate variation is as follows:
* for each lineage `l` in the network, the lineage's rate `r_l` is drawn from
  `lineage_distribution`, independently across lineages
Then, the length of each edge in each gene tree is multiplied by `r_l`.

The lineage rates `r_l` are saved as edge lengths in the network, in a new file
named `scalednetwork_lineagerates.phy` in the same directory as `genefile`
(see `save_lineagerates` below). `net` is *not* modified.

Output: gene trees, with edge lengths reflecting lineage rate variation,
(unless lineage_distribution = lognormal_meanone(0.0)),
without any degree-2 nodes unless `nodemapping=true,` written to `genefile`

# example

```repl
julia> net, level = simulatenetwork(4; seed=123, levelmax=1, verbose=false);

julia> writeTopology(net, round=true)
"(t2:0.448,((t1:0.202,#H10:0.0::0.06)I1:0.105,(t4:0.202,(t3:0.202)#H10:0.0::0.94)I2:0.105)I3:0.141)I4;"

julia> Random.seed!(520)

julia> gtrees = simulategenetrees_coalescentunits(net, 3, "trial_genetrees", lineage_distribution = lognormal_meanone(0.1));

julia> writeTopology(gtrees[3], round=true)
"((t1:0.88,t2:0.928):0.52,t3:1.398,t4:1.444);"

julia> Random.seed!(520)

julia> gtrees = simulategenetrees_coalescentunits(net, 3, "trial_genetrees"; lineage_distribution = lognormal_meanone(0.1), nodemapping=true);

julia> writeTopology(gtrees[3], round=true)
"(((t2:0.434)I4:0.273,((((t4:0.194)I2:0.09)I3:0.143)I4:0.056,(((t1:0.191)I1:0.106)I3:0.143)I4:0.056):0.217):0.031,((((t3:0.214)H10:0.0)I1:0.106)I3:0.143)I4:0.304);"

```
"""
function simulategenetrees_coalescentunits(net, nloci, genefile;
      lineage_distribution, nodemapping=false, nsites=1)
    gtcu = genefile * "_coal_unit" # before rate variation
    treelist = simulatecoalescent(net, nloci, 1; nodemapping=true)
    writeMultiTopology(treelist, gtcu) # gene trees with lengths in coalescent units
    length(treelist) == nloci || @warn "unexpected # of gene trees" # sanity check
    lineage_rate = Dict(e.number => rand(lineage_distribution) for e in net.edge)
    # add rate for population above the network's root. find its number first.
    rootedgenumber = max(0, maximum(e.number for e in net.edge)) + 1
    push!(lineage_rate, rootedgenumber => rand(lineage_distribution))
    open(genefile, "w") do fout
        for tree in treelist
            # add rate variation across lineages
            for e in tree.edge
                e.length *= lineage_rate[e.inCycle]
            end
            # cleanup gene trees
            # with PhyloNetworks v0.15: the degree-2 root is removed too, but that's fine for seq-gen
            nodemapping || PhyloNetworks.removedegree2nodes!(tree)
            write(fout, "[$nsites]")
            writeTopology(tree, fout)
            write(fout, "\n")
        end
    end  # close file cleanly
    netfile = joinpath(dirname(genefile), "scalednetwork_lineagerates.phy")
    save_lineagerates(net, lineage_rate, netfile) # saved inside edge lengths
    return treelist
end

"""
    function simulategenetrees(net::HybridNetwork, nloci::Int64, genefile::AbstractString, 
        nind::Int64, popsizes=nothing;
        lineage_distribution, nodemapping=false, nsites::Int64=1,
        substitutions_pergen::AbstractFloat = 1.25e-8,
        ingroupcrownheight = nothing)

Simulate `nloci` gene trees each with `nind` individuals per taxon along `net`
using PhyloCoalSimulations.
- if `popsizes` is nothing, then edge lengths are assumed to be in coalescent units
- if `popsizes` is a dictionary, then edge lengths are assumed to be in
  *number of generations*, and population sizes are used to simulate the coalescent.

The gene trees with branch lengths in the same units (before any scaling)
are written to file "`genefile`_networkunits".
These gene trees have extra degree-2 nodes to map them into the species network.

The model for rate variation is as follows:
* for each lineage `l` in the network, the lineage's rate `r_l` is drawn from
  `lineage_distribution`, independently across lineages
* all edges in all gene trees are then multiplied by the same scalar.

The lineage rates `r_l` are saved as edge lengths in the network, in a new file
named `scalednetwork_lineagerates.phy` in the same directory as `genefile`
(see `save_lineagerates` below). `net` is *not* modified.

Then, the length of edges in gene trees are scaled by an amount that
depends on the network branch that the gene tree maps to:
`r_l * scaler`
where `scaler = ngeneration_scaler * substitutions_pergen` to obtain
edge lengths in substitutions/site.

`ngeneration_scaler` is 1 (no effect) if `ingroupcrownheight` is nothing.
Otherwise, `ngeneration_scaler = ingroupcrownheight / ingrouptreeheight`
where `ingrouptreeheight` is calculated from the network as the distance
between the crown of the ingroup clade in the major tree and the leaves.
This option is meant to be used when the network is given with edge lengths
in coalescent units, and with an input `ingroupcrownheight` given
in number of generations, so that `ngeneration_scaler` scales the
coalescent units into a desired number of generations.
This implicitly determines the population size, to find the corresponding
number of generations needed to use the input rate in substitutions/generations.

Other output file:
gene trees written to `genefile_nsubst`, with edge lengths reflecting
lineage rate variation, in units of substitutions/sites
(no variation if lineage_distribution has 0 variance, e.g. lognormal_meanone(0.0)),
without any degree-2 nodes unless `nodemapping=true`.

```repl
julia> net, level = simulatenetwork(4; seed=123, levelmax=1, verbose=false);

julia> writeTopology(net, round=true)
"(t2:0.448,((t1:0.202,#H10:0.0::0.285)I1:0.105,(t4:0.202,(t3:0.202)#H10:0.0::0.715)I2:0.105)I3:0.141)I4;"

julia> Random.seed!(520)

julia> Ne = Dict(5  => 1022, 4  => 1366, 6  => 1061, 7  => 1446, 2  => 1152,
        10 => 1047, 9  => 1225, 8  => 1367, 3  => 1359, 1  => 1157)

julia> gtrees = simulategenetrees(net, 3, "trial_genetrees", 1, Ne, lineage_distribution = lognormal_meanone(0.0), substitutions_pergen = 1.25e-4);

julia> writeTopology(gtrees[3], round=true)
"((t1:0.006,t3:0.006):0.008,t2:0.014,t4:0.552);"

julia> gtrees = simulategenetrees(net, 3, "trial_genetrees", 2, Ne; lineage_distribution = lognormal_meanone(0.3), nodemapping=true, substitutions_pergen = 1.25e-4);

julia> writeTopology(gtrees[3], round=true)
"(((t2_2:0.009,t1_1:0.005):0.093,t3_1:0.093):0.45,((((t4_2:0.001,t1_2:0.001):0.004,t2_1:0.006):0.005,t4_1:0.007):0.007,t3_2:0.016):0.534);"
```
"""
function simulategenetrees(net::HybridNetwork, nloci::Int64, genefile::AbstractString, 
    nind::Int64, popsizes=nothing;
    lineage_distribution, nodemapping=false, nsites::Int64=1,
    substitutions_pergen::AbstractFloat = 1.25e-8,
    ingroupcrownheight = nothing)
  # file to save the gene trees with edge lengths in same units as in the network
  # (coalescent or # generations) before rate variation, and before x substitutions/site
  gtoriginal = genefile * "_networkunits" # 
  if isnothing(popsizes)
    treelist = simulatecoalescent(net, nloci, nind; nodemapping=true)
  else
    treelist = simulatecoalescent(net, nloci, nind, popsizes; nodemapping=false)
  end
  writeMultiTopology(treelist, gtoriginal) # gene trees with lengths in original units
  length(treelist) == nloci || @warn "unexpected # of gene trees" # sanity check
  # simulate a rate for each lineage in the network, including edge above root node
  # the input distribution should have mean 1: no units.
  lineage_rate = Dict(e.number => rand(lineage_distribution) for e in net.edge)
  rootedgenumber = max(0, maximum(e.number for e in net.edge)) + 1 # number of edge above root
  push!(lineage_rate, rootedgenumber => rand(lineage_distribution))
  # find height of the ingroup crown, which may be in coalescent units:
  # to later scale lengths to desired height, which may be in number of generations
  if isnothing(ingroupcrownheight)
    ngeneration_scaler = 1.0 # fixed network: keep the ingroup height (# generations) as is
  else
    #= adjust the overall edge lengths in case the network was randomly generated:
    if so, we want to control the total # of generations from the ingroup crown
    to the ingroup taxa to be the same across all generated networks.
    =#
    tree = majorTree(net)
    outgroups = outgrouplabels(net)
    for taxon in outgroups
        deleteleaf!(tree,taxon)
    end
    # now 'tree' is the major tree for the ingroup taxa
    ingrouptreeheight = maximum(PhyloNetworks.getHeights(tree))
    ngeneration_scaler = ingroupcrownheight / ingrouptreeheight
  end
  subrate_scaler = ngeneration_scaler * substitutions_pergen
  open(genefile* "_nsubst", "w") do fout
      for tree in treelist
          # add rate variation across lineages, and x substitutions/gen
          for e in tree.edge
              e.length *= lineage_rate[e.inCycle] * subrate_scaler
          end
          nodemapping || PhyloNetworks.removedegree2nodes!(tree)
          write(fout, "[$nsites]")
          writeTopology(tree, fout)
          write(fout, "\n")
      end
  end
  netfile = joinpath(dirname(genefile), "scalednetwork_lineagerates.phy")
  save_lineagerates(net, lineage_rate, netfile) # saved inside edge lengths
  return treelist
end

"""
    save_lineagerates(net, lineage_rate, netfile)

Save the rates in `lineage_rate` in file `netfile`, as edge lengths in `net`
in two ways, resulting in writing 2 networks:
1. edge length = lineage rate * original length in `net` (in coal units)
2. edge length = lineage rate itself
"""
function save_lineagerates(net, lineage_rate, netfile)
    netcp = deepcopy(net)
    for e in netcp.edge e.length *= lineage_rate[e.number]; end # rate * length
    writeTopology(netcp, netfile)
    for e in netcp.edge e.length  = lineage_rate[e.number]; end # rate itself only
    writeTopology(netcp, netfile; append=true)
    return nothing
end

"""
    lognormal_meanone(σ)

Lognormal distribution with mean one on the positive scale and
standard deviation σ on the log scale. The mean, on the log scale,
must be μ = -σ^2/2 for the mean to be 1 on the positive scale.
"""
function lognormal_meanone(sigma)
    mu = -sigma^2 / 2
    return LogNormal(mu, sigma)
end

"""
    function sim_snps(inputgts::AbstractString; nsites::Int64, ntries::Int64, 
        nbiallelic::Int64, juliaseed::Int64, 
        outputvcf::Abstract String, input_indls = ["t1", "t2", "t3", "t4"])

Given string path to a file of genetrees in Newick format, `inputgts`, and string vector
of individuals in each gene tree, `input_indls`, this function will simulate a SNP per 
gene tree (if the topology results in a sequence with variation).
For each gene tree, Seq-Gen will simulate a sequence of `nsites` length, taking the first
multiallelic site (2+ alleles). Seq-Gen simulates with set parameters:
    `-m HKY`: HKY model of nucleotide substitution- different rate of 
        transitions and transversions as well as unequal frequencies 
        of the four nucleotides (base frequencies)
    `-f 0.2 0.3 0.3 0.2`: base frequences of A, C, G and T respectively
    `-t 3`: transition/transversion ratio
    `-g 10`: the number of categories for the discrete gamma rate heterogeneity model
    `-a 0.35`: a real number >0 that specifies the shape of the gamma distribution to use with gamma rate heterogeneity
    `-of`: fasta output
    `-q`: quiet
    Parameters for Seq-Gen (-m, -a, -t, -g, -f) all adapted to reflect empirical data from the reptiles project.

If there are no multiallelic sites at first, Seq-Gen will simulate up 
to `ntries` times until a multiallelic (2+) site is produced, with a new seed per try produced 
from Julia seed `juliaseed`.  If after `ntries` times no multiallelic site is produced, 
that gene tree did not produce a SNP. `sim_snps` will stop simulating sequences with Seq-Gen
if requested number of biallelic SNPs `nbiallelic` is met. To determine whether `sim_snps`
met this requested number, check `seqgen_log`.

Output: VCF file with name/path `outputvcf`, with # biallelic SNPs not exceeding `nbiallelic`
(ideally would match `nbiallelic`, but depending on input gene trees, may not generate that many
SNPs). `outputvcf` will not have any monoallelic sites, but can have sites with 2+ alleles.
2 logs are also produced with information:
- `gene_log`: for genes that produced multiple alleles (2+), how many invariable sites needed 
        to be generated before the SNP was found.
- `seqgen_log`:contains information on how many genetrees were given as input (`input_genes`),
        how many biallelic SNPs were requested (`req_bi`), how many biallelic SNPs were successfully 
        simulated (`sim_bi`), how many multiallelic SNPs (3 or 4 alleles) were simulated (`sim_multi`,
        written to VCF, but don't count towards `sim_bi`), and how many genes failed to produce a SNP
        after `ntries` (`failed_genes`).
These will be outputted in the same directory as the `outputvcf`.
"""

function sim_snps(inputgts::AbstractString; nsites::Int64, ntries::Int64, 
        nbiallelic::Int64, juliaseed::Int64, 
        outputvcf::AbstractString, input_indls = ["t1", "t2", "t3", "t4"])
    blocklength = 10000
    #juliaseed = julia random number gen
    Random.seed!(juliaseed)
    #initialize counter of successful SNPs (some genes might not have any variation in sequence == unsuccessful)
    nsuccessfulsnp = 0
    counter_biallelic = 0
    #multi in this case == 3 or 4 alleles
    counter_multi = 0
    num_gts = countlines(inputgts)
    counter_failed = 0
    
    #check whether output vcf name contains a path
    #if so, put the seqgen logs in the same location
    if isdir(dirname(outputvcf))
        dir = dirname(outputvcf)
    #if not, use current directory as location for logs
    #this will be a problem if you are simulating multiple reps (logs will get overwritten)
    else dir = "." end

    #start logging
    logio = open("$dir/seqgen_log.txt", "w+")
    bygeneio = open("$dir/gene_log.txt", "w+")
    write(bygeneio, "gene_number,ninvar_sites\n")
    logger = SimpleLogger(logio)
    defaultlogger = global_logger(logger)
    
    open(inputgts, "r+") do trees_f
        position = -blocklength + 1
        genome_length = blocklength*nsites
        #get tab delimited string of indls 
        spstring = join(input_indls, "\t")
        #initialize pattern (minimizing garbage collecting)
        pattern = Dict{String, DNA}(sp => DNA_A for sp in input_indls)
        sequences = Dict{String,BioSequences.BioSequence}(sp => LongDNA{4}("") for sp in input_indls)
            
        open(outputvcf, "w+") do vcf_f
            #write header of vcf once, interpolating individuals and genome length
            write(vcf_f,"""
                ##fileformat=VCFv4.1
                ##contig=<ID=1,length=$genome_length>
                ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$spstring
                """)
            #make temporary directory within seqgen_temp dir
            isdir("seqgen_temp") || mkdir("seqgen_temp")
            dir = mktempdir("seqgen_temp", cleanup=true)
            
            #keeping track of gene tree number for log
            tree_num = 0
            #iterate over each gene tree
            for line in eachline(trees_f)
                tree_num += 1
                #write each indl gene tree to a file (necessary for seqgen input)
                #written in temporary directory
                open("$dir/gt", "w") do file write(file, line); end
                successfulgene = false
                #initialize tries to -1 -> first try will == 0
                    #then counter_invarsites can be j + tries*nsites
                tries = -1
                    
                for itry in 1:ntries
                    tries += 1   
                    sgseed = rand(1:10^20)
                    #run seq-gen with given nsites, seed, and input gene tree
                    sgcmd = `$seqgen -m HKY -f 0.2 0.3 0.3 0.2 -t 3 -g 10 -a 0.35 -l $nsites -of -z $sgseed -q $dir/gt`;
                    
                    #create io buffer for seq-gen output
                    seqgen_io = IOBuffer();

                    #if debugging seq-gen, change stderr to seqgen_io to write to log
                        #otherwise, lots of warnings clog up log when it's running fine
                        #(warning = "WARNING: The treefile contained partion lengths but only one partition was specified.")
                    run(pipeline(sgcmd, stdout = seqgen_io, stderr=devnull));
                    
                    seekstart(seqgen_io)
                    #create reader from seqgen output 
                    reader = FASTA.Reader(seqgen_io)
                    #collect records (identifier/sequence pairings) in vector (a record = an item)
                    record_vector = collect(reader)
                    #close reader
                    close(reader)
                    close(seqgen_io)

                    #create dictionary from records
                    sequences = Dict{String,BioSequences.BioSequence}(string(FASTA.identifier(record)) => LongDNA{4}(FASTA.sequence(record)) for record in record_vector)
                            
                    #get seq lengths, check if they are equal
                    all(length.(values(sequences)) .== nsites) || @error "sequences aren't of length $nsites\n"
                    
                    #check whether sequence has variable sites
                    if allequal(values(sequences)) 
                        continue #if no variable sites, resimulate
                    else #has variable sites, successfulgene== TRUE, break
                        successfulgene = true
                        nsuccessfulsnp += 1
                        break 
                    end
                end #ntries attempting to simulate sequences with variable sites

                if !successfulgene 
                    counter_failed += 1
                    continue
                end # to the next line, that is, next simulated gene tree
                    
                #if sequences have variable sites, loop until first sequence found
                j = 0 # index in the alignment
                #iterate over length of first sequence (having checked seqs all of same length)
                for _ in 1:nsites
                    j += 1
                    # (key, value)
                    for (sp, seq) in sequences
                        try 
                            pattern[sp] = seq[j]
                        catch e
                            @error "Got an empty sequence for a taxa/indl name, you probably didn't give the right input_indls\n"
                        end
                    end
                    #stop looping over each site once you find a site w/ variation
                    !allequal(values(pattern)) && break
                end
                counter_invarsites = j + tries*nsites
                #write # invarsites needed until a variable site is found to log
                write(bygeneio, "gene$tree_num,$counter_invarsites\n")

                #getting values/correctly formatted values for the vcf
                position += blocklength
                #get unique bases in pattern (max 4, min 2)
                uniquebases = unique(values(pattern))
                #increment biallelic counter if just 2 unique bases
                #increment multi counter if otherwise (should be 3 or 4 unique bases)
                length(uniquebases) == 2 ? (counter_biallelic += 1) : (counter_multi += 1)
                #reference base is first of unique bases
                refbase = string(uniquebases[1])
                #alt base becomes comma delimited string (if more than 2)
                altbase = join(uniquebases[2:end], ",")
                #initialize taxa_bases string (to be inputted to vcf) as empty
                taxa_bases = ""
                #create a dictionary matching unique bases to numeric values (starting at 0)
                base_dict = Dict(b => j-1 for (j, b) in enumerate(uniquebases))
                #get (string) of numeric values coding bases, tab delimited for each taxa
                for sp in input_indls
                    taxa_bases *= string(base_dict[pattern[sp]]) * "\t"
                end
                
                #write line per SNP (= per gene tree topology)
                write(vcf_f,"""
                        1	$position	.	$refbase	$altbase	.	.	.	GT	$taxa_bases
                        """)
                #fixit: check correct order of taxa 

                counter_biallelic >= nbiallelic && break #don't go to next gene tree if enough biallelic SNPs
                
            end #for each line in inputgts (each gene tree/ each resulting single SNP per GT)
            close(bygeneio)
        end #outputvcf do open/close block
        
        rm(dir, recursive=true) #delete temp dir because `cleanup=true` is doing it too late
    end #inputgts do open/close block
    
    #header for log
    write(logio, "input_genes,req_bi,sim_bi,sim_multi,failed_genes\n")
    #write stats from run
    write(logio, "$num_gts,$nbiallelic,$counter_biallelic,$counter_multi,$counter_failed\n")
    close(logio)
    global_logger(defaultlogger)
end

"""
    outgrouplabels(net::HybridNetwork, directedges=true::Bool;
                    subtreeroot = net.node[net.root],
                    ingroup_also = false::Bool)

Set of labels for tips in the outgroup, which is chosen to be the
descendants of the (first) edges of the root with the minimum number
of descendants.

# example
```julia
julia> net = readTopology("((hg19,ornAna1),(#H9:::0.4,((chrPic0,(allMis0,((galGal3,taeGut1))#H9:::0.6)),(anoCar2,pytMol0))));");

julia> outgrouplabels(net)
Set{String} with 2 elements:
  "ornAna1"
  "hg19"

julia> outgrouplabels(net; ingroup_also=true)
(["hg19", "ornAna1"], ["chrPic0", "allMis0", "galGal3", "taeGut1", "anoCar2", "pytMol0"])

julia> outgrouplabels(net; subtreeroot=net.node[16], ingroup_also=true) # outgroup subset of ingroup!
(["galGal3", "taeGut1"], ["chrPic0", "allMis0", "galGal3", "taeGut1", "anoCar2", "pytMol0"])

julia> outgrouplabels(net; subtreeroot=net.node[15], ingroup_also=true)
(["anoCar2", "pytMol0"], ["chrPic0", "allMis0", "galGal3", "taeGut1"])

julia> outgrouplabels(net; subtreeroot=net.node[9]) # this subtreeroot is a hybrid node
ERROR: the (subnetwork) root has a single daughter

julia> rootatnode!(net, "allMis0");

julia> outgrouplabels(net)
Set{String} with 1 element:
  "allMis0"

```
"""
function outgrouplabels(net::HybridNetwork, directedges=true::Bool;
                    subtreeroot::PhyloNetworks.Node = net.node[net.root],
                    ingroup_also::Bool = false)
    directedges && directEdges!(net)
    daughteredges = [e for e in subtreeroot.edge if PhyloNetworks.getParent(e)===subtreeroot]
    length(daughteredges) > 1 || error("the (subnetwork) root has a single daughter")
    if ingroup_also && length(daughteredges)>2
        @error("the (subnetwork) root has more than 2 daughters: no single ingroup. Turning off 'ingroup_only'.")
        ingroup_also = false
    end
    daughtersizes = [clustersize(ce) for ce in daughteredges]
    i = findmin(daughtersizes)[2] # first cluster of smallest size
    nodenum = PhyloNetworks.descendants(daughteredges[i])
    out_tiplabels = [n.name for n in net.node if n.number in nodenum]
    if ingroup_also
        i = 3-i # changes i: 2->1 or 1->2
        nodenum = PhyloNetworks.descendants(daughteredges[i])
        in_tiplabels = [n.name for n in net.node if n.number in nodenum]
        return (out_tiplabels, in_tiplabels)
    end
    return out_tiplabels
end

"""
    clustersize(edge::Edge)

Size of the edge's hardwired cluster (set of descendants).

`edge` should belong in a rooted network for which `isChild1` is up-to-date.
Run `directEdges!` beforehand on this network.
This is very important, otherwise one might enter an infinite loop,
and the function does not test for this.

# example

```julia
julia> net = readTopology("((hg19,ornAna1),(#H9:::0.4,((chrPic0,(allMis0,((galGal3,taeGut1))#H9:::0.6)),(anoCar2,pytMol0))));");

julia> clustersize(net.edge[3])
2

julia> clustersize(net.edge[17])
6

julia> # printEdges(net) # to see the list of edges with their parent & child nodes
```
"""
function clustersize(edge::PhyloNetworks.Edge)
    visited = Int[]
    ntips = clustersize!(visited, edge)
    return ntips
end
function clustersize!(visited::Vector{Int}, edge::PhyloNetworks.Edge)
    ntips = 0
    n = PhyloNetworks.getChild(edge)
    if n.hybrid # only need to check previous visits for hybrid nodes
        n.number in visited && return ntips
        push!(visited, n.number)
    end
    if n.leaf
        ntips += 1
    end
    for ce in n.edge
        if n == PhyloNetworks.getParent(ce)
            ntips += clustersize!(visited, ce)
        end
    end
    return ntips
end

"""
    function make_indfile(ind_names::Vector{String}, output_ind::AbstractString)

Make map file for individuals to populations, assuming individual names follow the convention
from `simulatecoalescent` (as implemented in `simulategenetrees`), where the original taxon name
is appended with the individual number (ex. if two individuals per taxa, "t1" -> "t1_1" and "t1_2").
Individuals are assumed to be in the same population if they share the same name before the underscore,
which becomes the population name.

Input: 
    - `ind_names`: vector of strings (individual names). Can be obtained
      with tipLabels(gt), where gt is a gene tree
    - `output_ind`: string name/location of output file.

Output: map file with name `output_ind`.
    Map file is created following the format:
    `individual U population`
    where U = unknown sex, and there is an individual per line.

#fixit: add note in docstring about underscore convention, if there are underscores previously, they are not a problem (confirm)    
"""

function make_indfile(ind_names::Vector{String}, output_ind::AbstractString, re= r"_\d+$")
#fixit just double check this works
    species = Dict()
    for ind in ind_names
        pop = replace(ind, re => "")
        species[ind] = pop   
    end

    open(output_ind, "w+") do io
        for (key, value) in species
            write(io, "$key\tU\t$value\n")
        end
    end
end

"""
    parameterset_directory(net, nind, ngenes,
        nbiallelic, lineage_dist, sub_rate)

Directory for output files for all replicates given a set of parameters.
Used by several functions, defined once here for consistency across the pipeline.
"""
function parameterset_directory(net, nind, ngenes,
        nbiallelic, lineage_dist, sub_rate)
    joinpath("output",
      join(["$net", "$nind-indls", "$ngenes-genes", 
            "$nbiallelic-biall", "$lineage_dist-lindist", "$sub_rate-subrate"],
            "_")
    )
end

"""
    readmultigts_seqgen(file::AbstractString)

Modified version of `readMultiTopology` for use on gene tree files written in 
seq-gen format, ex.
```
[100](((t2_2:0.52,t4_1:2.42):0.55,t2_1:0.67):0.84,((t1_1:0.47,t4_2:2.46):1.95,t3_2:0.85):0.69,(t3_1:0.23,t1_2:0.68):0.43);
```
with the number of loci in brackets before the Newick formatted gene tree.

This function strips the loci number for each line, and reads the Newick gene tree into a vector of `HybridNetwork`s.
"""

function readmultigts_seqgen(file::AbstractString)
    s = open(file)
    numl = 1
    vnet = HybridNetwork[];
    for line in eachline(s)
        line = strip(line) # remove spaces
        line = split(line, "]")[2]
        c = isempty(line) ? "" : line[1]
        if(c == '(') push!(vnet, readTopology(line,false)) end # false for non-verbose
        numl += 1
    end
    close(s)
    return vnet
end

"""
    simulate_fromnet_tovcf(net::HybridNetwork; nind::Int, 
        lineage_dist::AbstractFloat, sub_rate::AbstractFloat, popsizes=nothing,
        nsites::Int, nbiallelic::Int, ngenes::Int=round(Int, 1.5*nbiallelic), ntries::Int=20, 
        nrep::Int, rep_id_start::Int, seed::Int=nothing, netname="fleg")

Given an input network `net`, this function:
1. simulates `ngenes` (default = 1.5*nbiallelic) number of gene trees within that network,
   with `nind` number of individuals simulated per taxon and substitutions/gen `sub_rate`. 
   Lineages across gene trees are scaled pulling rates from the lognormal 
   distribution with mean one and standard deviation `lineage_dist`.
   Output: `gts_nsubst`, `gts_networkunits`, `scalednetwork_lineagerates.phy`
2. simulates SNPs from gene trees, by simulating a sequence of `nsites` length up to
   `ntries` times (default=20). The first multiallelic site (2+ alleles) from the sequence will be saved to
   a VCF. Once `nbiallelic` number of biallelic SNPs are simulated, `sim_snps` will stop
   simulating a sequence per gene tree. Check logs to see how many biallelic vs. 
   multiallelic (3+) SNPs were generated, as well as how many invariable sites it took per
   gene to generate a SNP.
   Output: `seq.vcf`, `gene_log.txt`, `seqgen_log.txt` 
   `gene_log.txt` has a line for each gene/SNP simulated, with the number of invariable sites simulated
        until a variable site was found for that gene.
    `seqgen_log.txt` contains information on how many genetrees were given as input (`input_genes`),
        how many biallelic SNPs were requested (`req_bi`), how many biallelic SNPs were successfully 
        simulated (`sim_bi`), how many multiallelic SNPs (3 or 4 alleles) were simulated (`sim_multi`,
        written to VCF, but don't count towards `sim_bi`), and how many genes failed to produce a SNP
        after `ntries` (`failed_genes`).
3. creates a file mapping individuals to populations, based on the naming scheme in
   PhyloCoalSimulations (ex. appending a `_1` or `_2` to the taxon name if there are
   two individuals). The file may have individuals listed in a different order than `seq.vcf`,
   but this is not a problem, as individual ordering in the indmap doesn't matter for VCF 
   to eigenstrat converstion.
   Output: `indmap`
4. converts the VCF file to eigenstrat format, using the `indmap`. Order of individuals will match the
   input VCF file, not the input `indmap` file.
   Output: `seq.geno`, `seq.ind`, `seq.snp`

This function is intended for simulating in replicate. `nreps` number of replicates will be 
simulated, starting from `rep_id_start` starting number. Optional argument `seed`

See example use in file `simulate_batch.jl`, where replicates are parallelized.
"""

function simulate_fromnet_tovcf(net::HybridNetwork; nind::Int, 
        lineage_dist::AbstractFloat, sub_rate::AbstractFloat, popsizes=nothing,
        nsites::Int, nbiallelic::Int, ngenes::Int=round(Int, 1.5*nbiallelic), ntries::Int=20, 
        nrep::Int, rep_id_start::Int, seed::Int=nothing, netname="fleg")

    #get starting datetime
    datetime = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

    net_str = netname
    randomnet = (netname != "fleg")
    #boolean that's false if netname is "fleg"
        #true if netname is not "fleg"

    #create folder for output if it doesn't exist already
    paramset_dir = parameterset_directory(net_str, nind, ngenes,
        nbiallelic, lineage_dist, sub_rate)
    isdir(paramset_dir) || mkdir(paramset_dir)
    isnothing(seed) || Random.seed!(seed)
    # simulate all rep seeds, in case resps need to be re-run in different order
    seqgen_seed_vector = rand(1:10000, rep_id_start + nrep)

    rep_id = rep_id_start
    rep_id_end = rep_id_start + nrep
    while rep_id < rep_id_end
        # create subfolder for this rep if it doesn't exist already
        rep_dir = joinpath(paramset_dir, string(rep_id, pad = 3))
        isdir(rep_dir) || mkdir(rep_dir)
            
        #now with successful rep, start logging
        io = open("$rep_dir/log.txt", "a+")
        logger = SimpleLogger(io)
        defaultlogger = global_logger(logger)
        write(io, """
        $datetime
        seed = $seed
        """)

        #simulate gts, write to rep_dir folder
        if isfile("$rep_dir/gts_nsubst")
            gts = readmultigts_seqgen("$rep_dir/gts_nsubst")
            ind_names = tipLabels(gts[1])
            println(io, "gene trees already exist")
        else
            gts = simulategenetrees(net, ngenes, "$rep_dir/gts", nind, popsizes;
                lineage_distribution = lognormal_meanone(lineage_dist), 
                nodemapping=false, nsites=nsites,
                substitutions_pergen = sub_rate,
                ingroupcrownheight = nothing)
            ind_names = tipLabels(gts[1])
            write(io, """
            nind = $nind
            ngenes = $ngenes
            lineage_dist = $lineage_dist
            sub_rate = $sub_rate
            nsites = $nsites
            nbiallelic = $nbiallelic
            """)
        end

        # simulate sequences with seq-gen
        if isfile("$rep_dir/seq.vcf")
            println(io, "sequence already exists")
        else
            sequence_seed = seqgen_seed_vector[rep_id]
            @show(sequence_seed)
            println(io, "sequence_seed = $sequence_seed")
            sim_snps("$rep_dir/gts_nsubst"; nsites=nsites, ntries=ntries, 
                nbiallelic=nbiallelic, juliaseed=sequence_seed, 
                outputvcf= "$rep_dir/seq.vcf", input_indls = ind_names)
        end
        
        if isfile("$rep_dir/indmap")
            println(io, "indmap already exists")
        elseif nind == 1
            #just map each indl -> pop
            #pass re as "" (empty string) to make_indfile
            make_indfile(ind_names, "$rep_dir/indmap", "")
        else #if nind > 1, map individual names to populations and make indmap file
            make_indfile(ind_names, "$rep_dir/indmap")
        end

        if isfile("$rep_dir/seq.geno") && isfile("$rep_dir/seq.snp")
            println(io, "eig files already exist")
        else
            run(`python3 scripts/vcf2eigenstrat.py -v $rep_dir/seq.vcf -o $rep_dir/seq -i $rep_dir/indmap`)
        end

        @show(rep_id)
        rep_id += 1
        close(io)
        global_logger(defaultlogger)

    end
end

println("done reading this awesome file")