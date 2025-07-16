#= code to better understand the similarities and differences of the 14 graphs
at hardwired-distance 0 from G4, yet non-isomorphic to G4.

see `summarize_graphs.qmd` subsection "Debugging hash matches" for how these
graphs were identified.
=#
using CSV
using DataFrames
using TidierData
using PhyloNetworks
using PhyloPlots
using RCall

# load fleg (with bottleneck node) and fleg_major (bottleneck node suppressed)
include("input/fleg-net.jl")
nearmajortree = displayedtrees(fleg, 0.45) # still has bottleneck node
removedegree2nodes!.(nearmajortree)  # removes the root also

alldisplayed = removedegree2nodes!.(displayedtrees(fleg, 0.0))
# 16 trees: but some have the same topology
for j in reverse(1:16)
  for i in 1:(j-1)
    if hardwiredclusterdistance(alldisplayed[j], alldisplayed[i], true) == 0
      deleteat!(alldisplayed, j)
      break # don't try other i's!
    end
  end
end
length(alldisplayed) # 8 displayed tree *topologies* only

function ndisplayed_intargetlist(net::HybridNetwork, trees2find)
  dtrees = removedegree2nodes!.(displayedtrees(net, 0.0))
  ndisplayed_intargetlist(dtrees, trees2find)
end
function ndisplayed_intargetlist(candidatetrees::Vector{HybridNetwork}, trees2find)
  ntargets = length(trees2find)
  notfound = collect(1:ntargets)
  for tree in candidatetrees
    isempty(notfound) && break
    for (ii,i) in enumerate(notfound)
      if hardwiredclusterdistance(tree, trees2find[i], false) == 0
        deleteat!(notfound, ii) # no need to look for trees2find[i] in later displayed trees
        break # 'tree' cannot match other near-major trees
      end
    end
  end
  return ntargets - length(notfound)
end

# read the input graphs: at hardwired-distance 0 from the true G4
df_52 = CSV.read("output/g4_hwmatches.csv", DataFrame)
df_52 = @chain df_52 begin
    @mutate(n = row_number())
end
nets_52 = readMultiTopology(df_52[!,:newick])

# graph 14 on line 15 seems to be equal to G4:
# in the output file in which graph 14 was found, the network with topology G4
# and edge parameters optimized by qpGraph had this hash and newick string:
flegnobott_hash = "b66401d14d13e41dfd92363767147547"
# G4 with parameters optimized by qpGraph, and data for graph 14
"(chimp:0.0997828322130873,(((denisovan:0.014679814613136222)#H1:0.0::0.7451647456609543,((#H1:0.0::0.25483525433904575,(neanderthal_west:0.000143937289744886,(neanderthal_east:0.0019323896297988904,#H2:0.0::0.046915160325463286):0.0):0.01378273003610523):0.01379603107975491,((#H3:0.0::0.3739290032754603)#H2:0.0::0.9530848396745367,(((africa_east:0.001172536939233813,((non_africa_west:0.0008552315577179843,non_africa_east:0.0024381932503937127):0.006916024063897025)#H3:0.0::0.6260709967245397):0.003850218601766436,africa_west:0.0029024842825797347):0.0008445226429162036)#H4:0.0::1.0):0.001107295577989577):0.0008634867031558632):0.09978283221308759,#H4:0.0::0.0):0.0997828322130861);"

nets_g14hash = @chain df_52 begin
    @filter(hash == "b66401d14d13e41dfd92363767147547")
end
show(
  select(nets_g14hash, [:n,:hash,:truenet,:net_hwdist,:tree_hwdist,:ndisplayed,:irep,:irun]),
  allcols=true)
#=
5×8 DataFrame
 Row │ n      hash                              truenet  net_hwdist  tree_hwdist  ndisplayed  irep   irun  
     │ Int64  String                            Bool     Int64       Int64        Int64       Int64  Int64 
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │    14  b66401d14d13e41dfd92363767147547    false           0            0           2     20     32
   2 │    31  b66401d14d13e41dfd92363767147547    false           0            0           2     13     18
   3 │    35  b66401d14d13e41dfd92363767147547    false           0            2           2     51     21
   4 │    42  b66401d14d13e41dfd92363767147547    false           0            0           2     92      1
   5 │    48  b66401d14d13e41dfd92363767147547    false           0            0           2      1     16
Networks 14 and 35 have the same hash,
yet different tree_hwdist: they are isomorphic *if* we ignore the major/minor
attribute of hybrid edges.
Network 35 has the gene flow neanderthal stem -> denisovan as major instead of minor.
=#

#= Look at networks with same hash as G4 manually,
to decide if they are truly isomorphic to g4.
Conclusion: yes. For network 35, we need to ignore the major/minor attributes
=#
df_g4hash = subset(df_52, :hash => x -> x .== flegnobott_hash)
R"cairo_pdf"("tmp/fig_graphs_samehashas_g4.pdf", height=10, width=10)
R"par"(mar=[0,0,0,0]); R"layout"([1 2; 3 4; 5 6]);
for n in df_g4hash.n
  plot(nets_52[n], xlim=[1,13])
  R"mtext"("graph $n", line=-2)
end
R"dev.off"()

# Are there other examples of graphs with the same hash yet distinct tree_hwdist?
# yes: at least 5 examples from our file
unique(df_52[!,:net_hwdist]) # 0: that was the criterion to select these graphs
unique(df_52[!,:ndisplayed]) # 2: all graphs display the 2 nearmajor trees (?)
rangeof(x) = (mM = extrema(x); return mM[2] - mM[1])
tmp = @chain df_52 begin
  @group_by(hash)
  @summarize(
    n = first(n),
    numnets = n(),
    tree_hwdist_min = minimum(tree_hwdist),
    tree_hwdist_max = maximum(tree_hwdist))
  @filter(tree_hwdist_min < tree_hwdist_max)
end
#=
 Row │ hash                              n      numnets  tree_hwdist_min  tree_hwdist_max 
     │ String                            Int64  Int64    Int64            Int64           
─────┼────────────────────────────────────────────────────────────────────────────────────
   1 │ b6846e5befa755dba55317b77eeb409d      1        9                0                2
   2 │ 40e4366dd88fff3c6924481f8e2fdfb0      5        7                0                2
   3 │ fbdd199eafbcfcee4cc451783c727dbf     11       11                0                2
   4 │ b66401d14d13e41dfd92363767147547     14        5                0                2
   5 │ 38118a62079e2df765d364ece9187923     18        2                0                2

Each hash above is the hash of at least 1 graph with tree_hwdist 0,
and of a different graph with tree_hwdist 2.
We already looked at the 4th hash on this list:
graph 14 != graph 35 due to minor/major switch, yet same hash.
Here, all these graphs have hardwired cluster dissimilarity 0 with G4 (=graph 14).

To get a list subset:
=#

@chain subset(df_52, :hash => ByRow(x -> x in tmp.hash)) begin
  @group_by(hash, tree_hwdist)
  @summarize(n = first(n), numnets = n())
  @ungroup
  @arrange(hash, tree_hwdist)
end
#=
10×4 DataFrame
 Row │ hash                              tree_hwdist  n      numnets 
     │ String                            Int64        Int64  Int64   
─────┼───────────────────────────────────────────────────────────────
   1 │ 38118a62079e2df765d364ece9187923            0     27        1
   2 │ 38118a62079e2df765d364ece9187923            2     18        1
   3 │ 40e4366dd88fff3c6924481f8e2fdfb0            0      5        5
   4 │ 40e4366dd88fff3c6924481f8e2fdfb0            2     30        2
   5 │ b66401d14d13e41dfd92363767147547            0     14        4
   6 │ b66401d14d13e41dfd92363767147547            2     35        1
   7 │ b6846e5befa755dba55317b77eeb409d            0      1        2
   8 │ b6846e5befa755dba55317b77eeb409d            2      4        7
   9 │ fbdd199eafbcfcee4cc451783c727dbf            0     17        1
  10 │ fbdd199eafbcfcee4cc451783c727dbf            2     11       10

For example, the following pairs have the same hash, probably same "raw" topologies
but different minor/major hybrid edge attributes somewhere:
- graphs 18 & 27
- graphs  5 & 30
- graphs 14 & 35 (already seen)
- graphs  1 &  4
- graphs 11 & 17
=#

# keep first network from each hash
df_14 = @chain df_52 begin
    @group_by(hash)
    @summarise(newick = first(newick), numnets = n(), n = first(n))  
end
nets = readMultiTopology(df_14[!,:newick])

nn = nrow(df_14)
df_14.i = collect(1:nn)
df_14.truenet_hwd   = Vector{Union{Missing, Int}}(missing, nn)
df_14.majortree_hwd = Vector{Union{Missing, Int}}(missing, nn)
df_14.nnearmajor_displayed = Vector{Union{Missing, Int}}(missing, nn)
df_14.ndisplayed_displayed = Vector{Union{Missing, Int}}(missing, nn)
for (i,net) in enumerate(nets)
  df_14.truenet_hwd[i] = hardwiredclusterdistance(net, fleg_nobott, true)
  df_14.truenet_hwd[i] == 0 || @error("graph $i was supposed to be at distance 0...")
  df_14.majortree_hwd[i] = hardwiredclusterdistance(majortree(net), nearmajortree[1], false)
  dtrees = removedegree2nodes!.(displayedtrees(net, 0.0))
  df_14.nnearmajor_displayed[i] = ndisplayed_intargetlist(dtrees, nearmajortree)
  df_14.ndisplayed_displayed[i] = ndisplayed_intargetlist(dtrees, alldisplayed)
end

@chain df_14 begin
    @filter(majortree_hwd==0)
    @select(n, i, # hash, newick,
      numnets, truenet_hwd, majortree_hwd,
      nnearmajor_displayed, ndisplayed_displayed)
end

#= network 9 (after filtering, originally 14) is most similar:
- same major tree as in G4
- displays both nearmajor trees
- displays all 8 trees displayed by G4

Next most similar network: 4 after filtering, originally 5:
- same major tree as in G4
- displays both nearmajor trees
- displays 7 of the 8 trees displayed by G4
=#

# now look for least similar graph
minimum(df_14[!,:ndisplayed_displayed]) # 5
@chain df_14 begin
    @filter(ndisplayed_displayed==5)
    @select(n, i, # hash, newick,
      numnets, truenet_hwd, majortree_hwd,
      nnearmajor_displayed, ndisplayed_displayed)
end
# example: i=3 after filtering (originally 3 also)


# trial plots

mynet = deepcopy(nets[9]); # same as G4!!
for i in [-3,-4,-5,-6,-9,-15] rotate!(mynet, i); end

mynet4 = deepcopy(nets[4]);
for i in [-3,-5,-6,-8,-9,-10,-12] rotate!(mynet4, i); end
setgamma!(mynet4.edge[2], 0.51)  # was γ=0.31866503546018676

mynet3 = deepcopy(nets[3]);
for i in [-3,-4,-6,-7,-13, -9,-10] rotate!(mynet3, i); end

plot(mynet; showgamma=true, shownodenumber=true);
#= conclusion:
- nets[9] is isomorphic to G4
- nets[4] differs by the ghost gene flow into humans being placed
  after the split denisovan - neandethals + humans, instead of before.
- nets[3] differs by the donor of the later ghost into nonafricans:
  starting from the first ghost lineage instead of
  starting from the human stem (before receiving first ghost)
=#

# figure: similar settings (height, etc.) as for the G1 & G4 figures
#         in figures/intro_methods/fig-flegsidebyside.jl
# first: capitalize Neanderthal, Africa, Denisovan; remove underscores
for net in (mynet4, mynet3), n in net.leaf
  n.name = replace(n.name, "non_" => "non-", "_" => " ","chimp" => "chimpanzee",
    "nean" => "Nean", "africa" => "Africa", "deni" => "Deni")
end
R"cairo_pdf"("tmp/fig_graphs_dist0from_g4.pdf", height=5, width=6.5)
R"layout"([1 1; 2 2]);
R"par"(mar=[0,0,0,0], oma=[0,0,0,0]);
# 5th graph in g4_hwmatches.csv. has 7 of the 8 trees displayed by g4
plot(mynet4, xlim=[1,13], tipoffset=0.2, edgelabelcolor="red4", edgecex=0.6);
R"mtext"("a)", line=-2, adj=0.01);
# 3rd graph in g4_hwmatches.csv. has 5 of the 8 trees displayed by g4
plot(mynet3, xlim=[1,13], tipoffset=0.2, edgelabelcolor="red4", edgecex=0.6);
R"mtext"("b)", line=-2,adj=0.01);
# R"mtext"("Example networks with identical hardwired clusters as G4, yet distinct from G4",
#  outer=true, side=3, line=-1, cex=0.9)
R"dev.off"()
