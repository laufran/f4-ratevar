using PhyloNetworks
using DataFrames
using CSV
using Combinatorics

#=
miscellaneous functions, see their docstrings for details:

    nodenumber2taxonname(net, number)
    calc_quartet_numhybrids(net, csv_filename)
    net_to_adjlist(net)
    edgelist_to_net()

and helper functions without docstrings:

    quarnet_type(net, fourtaxa_vector)
    numhybrids_in_4blob(quarnet)
    quartetsplit(net, descendant_numbers, fourtaxonmames_vector)
=#

"""
    nodenumber2taxonname(net, num)

`node.name` for the node in `net` such that `node.number` is `num`.
"""
function nodenumber2taxonname(net, num)
    i = findfirst(n -> n.number == num, net.node)
    net.node[i].name
end

"""
    calc_quartet_numhybrids(net::HybridNetwork, outputcsv::AbstractString)

Given an input HybridNetwork `net`, for every four-taxon/individual subnetwork
possible in the network, calculate:
1. taxonset: the 4 taxa names sorted alphabetically, then concatenated into a
   string with individual names separated by dashes. example: `t1-t2-t3-t4`
2. split: empty string if the 4-taxon subnetwork is not "tree-like".
   The network is tree-like if it has a cut-edge that separates 2 taxa from
   the other 2 (these 2 pairs are disconnected when theh edge is cut).
   In this case, the split pair is a string of the form `t1-ta|tb-tc`
   (continuing on the previous example) where `ta` is either `t2`, `t3` or `t4`,
   and `tb < tc` are the 2 taxa other than `t1` and `ta`.
3. split code: "12_34", "13_24", or "14_23" if the 4-taxon subnetwork is tree-like,
   depending on which taxon pairs up with the first taxon in the split.
   Nothing if the 4-taxon subnetwork is not tree-like (does not have a 2-2 split).
4. h1 (integer): 0 if the 4-taxon subnetwork is tree-like. Otherwise, it
   has a "4-blob", that is, a biconnected component incident to 4 cut-edges,
   one leading to each of the 4 taxa. In that case, h1 is the number of
   hybridization events in that 4-blob,
   *without* shrinking 2-cycles
5. h2 (integer): 0 if the 4-taxon subnetwork is tree-like, and otherwise
   the number of hybridization events in the subnetwork's 4-blob
   after shrinking all 3-cycles
6. the gammas: at most 3 minor gamma value(s) in the 4-blob after shrinking 3-cycles,
   sorted in descending order. This will be empty if the 4-taxon network
   is tree-like (no 4-blob: h1=h2=0).

Each quartet will form a row in output csv `outputcsv` (string that forms the 
file name/output path). example format:

    taxonset,splitpair,h1,h2,gamma1,gamma2,gamma3
    t1-t3-t5-t6,t1-t5|t3-t6,13_24,2,1,0.279,NA,NA

# example

```julia
julia> include("input/fleg-net.jl"); # defines network "fleg"

julia> # R"par"(mar=[0,0,0,2]); plot(fleg, showgamma=true);

julia> include("scripts/miscfxns.jl")

julia> fourtaxa_df = calc_quartet_numhybrids(fleg, "fleg-quarnets.csv");

julia> fourtaxa_df[[1,4,8],:] # rows 1 (tree-like), 4 and 8 (blob-like)
3×8 DataFrame
 Row │ fourtaxa                           split                              splitcode  h1      h2      g1          g2          g3         
     │ String?                            String?                            String?    Int64?  Int64?  Float64?    Float64?    Float64?   
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ africa_east-africa_west-chimp-de…  africa_east-africa_west|chimp-de…  12_34           0       0  missing     missing     missing    
   2 │ africa_east-africa_west-chimp-no…                                                     3       1        0.32  missing     missing    
   3 │ africa_east-africa_west-denisova…                                                     4       4        0.52        0.32        0.19

julia> fourtaxa_df[[1,4,8], 1] # to see full taxon names
3-element Vector{Union{Missing, String}}:
 "africa_east-africa_west-chimp-denisovan"
 "africa_east-africa_west-chimp-non_africa_east"
 "africa_east-africa_west-denisovan-non_africa_east"

```
"""
function calc_quartet_numhybrids(net::HybridNetwork, outputcsv::AbstractString)
    taxonlist = sort(tiplabels(net))
    ntaxa = length(taxonlist)
    ntaxa >= 4 || error("fewer than 4 taxa, 0 quarnet")
    nquartets = binomial(ntaxa, 4)
    df = DataFrame(
        fourtaxa  = Vector{Union{Missing, String}}(missing, nquartets),
        split     = Vector{Union{Missing, String}}(missing, nquartets),
        splitcode = Vector{Union{Missing, String}}(missing, nquartets),
        h1 = Vector{Union{Missing, Int}}(missing, nquartets),
        h2 = Vector{Union{Missing, Int}}(missing, nquartets),
        g1 = Vector{Union{Missing, Float64}}(missing, nquartets),
        g2 = Vector{Union{Missing, Float64}}(missing, nquartets),
        g3 = Vector{Union{Missing, Float64}}(missing, nquartets),
    )
    for (j,fourtaxa) in enumerate(combinations(taxonlist, 4))
        df[j,:fourtaxa] = join(fourtaxa, '-')
        sp,h1,h2,γ = quarnet_type(net, fourtaxa)
        df[j,:split]     = sp[1]
        df[j,:splitcode] = sp[2]
        df[j,:h1] = h1
        df[j,:h2] = h2
        df[j,:g1] = γ[1]
        df[j,:g2] = γ[2]
        df[j,:g3] = γ[3]
    end
    CSV.write(outputcsv, df)
    return df
end

function quarnet_type(net::HybridNetwork, fourtaxa)
    length(fourtaxa) == 4 || error("there aren't 4 tips: $fourtaxa")
    # check that all 4 taxa are in net
    fourtaxa ⊆ tiplabels(net) ||  error("tips are not all in network")
    # create 4-taxon subnetwork
    quarnet = deepcopy(net)
    for tip in tiplabels(net)
        tip ∈ fourtaxa && continue
        deleteleaf!(quarnet, tip, simplify=false)
    end
    # make it semidirected: remove above LSA and unroot
    deleteaboveLSA!(quarnet)
    if length(getroot(quarnet).edge) == 2
      PhyloNetworks.fuseedgesat!(quarnet.rooti, quarnet)
    end
    removedegree2nodes!(quarnet) # flegontov has 1 degree-2 node for bottleneck
    # get h in blob before shrinking 3-cycles
    h1, split_descendant_numbers, gammas = numhybrids_in_4blob(quarnet)
    if h1 > 0
      shrink3cycles!(quarnet, true) # 'true' to unroot the network. also shrinks 2-cycles.
      h2, _, gammas = numhybrids_in_4blob(quarnet)
    else
        h2=0 # gammas unchanged: should all be missing
    end
    h1 ≥ h2 || error("""h1=$h1 < h2=$h2, 4 taxa: $fourtaxa
        quarnet: $(writenewick(quarnet, round=true))""")
    if h2 == 0
        all(ismissing.(gammas)) || error("h2=0 yet non-missing gammas")
        h1 == 0 || error("""
          h2=0 yet h1=$h1: the 4-blob disappeared after shrinking 2- and 3-cycles. wth?
          4 taxa: $fourtaxa
          quarnet: $(writenewick(quarnet, round=true))
          """)
    end
    # if h1=h2=0, quarnet is tree-like: find its split
    if isnothing(split_descendant_numbers)
        split = ["",""]
    else # find split expression: 12_34 etc.
        split = quartetsplit(quarnet, split_descendant_numbers, fourtaxa)
    end
    return split, h1, h2, gammas
end

# network `quarnet` should have exactly 4 taxa, not checked
function numhybrids_in_4blob(quarnet)
    blobs = biconnectedcomponents(quarnet) # default: includes trivial blobs
    articulations = Set{Int}() # numbers of a blob's articulation nodes. re-used
    # a 4-blob has 4 articulation nodes, if the network is binary
    for blobedges in blobs # 1 blob = 1 vector of edges
        empty!(articulations)
        length(blobedges) > 1 || continue # ignore trivial blobs
        for e in blobedges
            for node in e.node
                node.number in articulations && continue
                length(node.edge) == 3 || error("non-binary net at node $node.\nnet: $(writenewick(quarnet, round=true))")
                if !(node.edge ⊆ blobedges)
                    push!(articulations, node.number)
                end
            end
        end
        if length(articulations) == 4 # found 4-blob, which much be unique
            gammas = sort!([e.gamma for e in blobedges if !e.ismajor]; rev= true)
            h = length(gammas)
            largest3gammas = Vector{Union{Missing, Float64}}(missing, 3)
            for i in eachindex(gammas)
                largest3gammas[i] = gammas[i]
                i == 3 && break
            end
            return (h, nothing, largest3gammas)
        end
    end
    # no 4-blob was found: h = 0 and the quarnet is tree-like. What split?
    # find a cut edge with 2 descendants, return these descendants' node numbers
    emptygammas = Vector{Union{Missing, Float64}}(missing, 3)
    for blobedges in blobs
        length(blobedges) == 1 || continue # cut-edges only
        descendant_nodenumbers = PhyloNetworks.descendants(blobedges[1])
        number_descendants = length(descendant_nodenumbers)
        if number_descendants == 2
            return (0, descendant_nodenumbers, emptygammas)
        end
    end
    @error("""no 4-blob (h=0) yet no 2-2 split.
          quarnet: $(writenewick(quarnet, round=true))""")
    return (0, nothing, emptygammas)
end

# tax = vector of 4 taxon names, not checked
function quartetsplit(net, descendant_numbers, tax)
    descendant_names = [nodenumber2taxonname(net, num) for num in descendant_numbers]
    sbv = BitArray(n in descendant_names for n in tax) # split bit vector
    sum(sbv) == 2 || error("split is not 2-2, 4 taxa: $tax, descendants: $descendantnames")
    # make it one of 3, starting with 1: 1100, 1010, 1001
    if !sbv[1]
        sbv .= .!sbv
    end
    split_code = (sbv[2] ? "12_34" : (sbv[3] ? "13_24" : "14_23"))
    split = (sbv[2] ? "$(tax[1])-$(tax[2])|$(tax[3])-$(tax[4])" :
            (sbv[3] ? "$(tax[1])-$(tax[3])|$(tax[2])-$(tax[4])" :
                      "$(tax[1])-$(tax[4])|$(tax[2])-$(tax[3])"))
    return split, split_code
end

# fixit: consider using the version of these conversion functions from here:
# https://github.com/JuliaPhylo/PhyloUtilities/blob/admix/scripts/interop_admixtools.jl
"""
  net_to_adjlist(net::HybridNetwork)
  net_to_adjlist(net_newickfilename::AbstractString)

Given a `net` of either HybridNetwork format, string path to text file location,
or string Newick format, `net_to_adjlist` returns a two column DataFrame
with source/target node pairings.

```repl
julia> net = readnewick("((t4:1.504,((t3:1.268,t1:1.268):1.057,#H8:0.0::1.424):1.179):1.16,(t2:1.325)#H8:1.34::1.576);")
julia> adjlist = net_to_adjlist(net)
source target
1     i1     t3
2     i1     t1
3     i2     i1
4     i2     H8
5     i3     t4
6     i3     i2
7     H8     t2
8     i4     i3
9     i4     H8
```
"""
function net_to_adjlist(net::HybridNetwork)
    adj_list = DataFrame(source = String[], target= String[]) # initialized as empty
    PhyloNetworks.resetedgenumbers!(net)
    PhyloNetworks.resetnodenumbers!(net)
    PhyloNetworks.nameinternalnodes!(net, "i") # "i" = prefix for new names
    for n in net.node
        snode_name = n.name # source node name
        for c in getchildren(n) # will do nothing if leaf
            tnode_name = c.name # name of child (target) node
            push!(adj_list, [snode_name,tnode_name]) # push new row
        end
    end
    return adj_list
end

net_to_adjlist(netpath::AbstractString) = net_to_adjlist(readnewick(netpath))

"""
    edgelist_to_net(
        edgesource::AbstractVector,
        edgetarget::AbstractVector,
        edgetype::AbstractVector,
        edgeweight=missing,
        weighttol = -1e-8
    )

HybridNetwork object built from a list of edges, where each edge corresponds to a row.
The source is the name of the parent node;
the target is the name of the child node;
the type is the string "edge" for tree edges, "admix" for hybrid edges.
The weight is the edge length for tree edges and the inheritance probability for
hybrid edges. If not given, then all of them are given the default value of 0.5.

# example

```repl
julia> net0 = readnewick("((t4:1.504,((t3:1.268,t1:1.268):1.057,#H8:0.0::0.475):1.179):1.16,(t2:1.325)#H8:1.34::0.525);");

julia> al = net_to_adjlist(net0); esource = al[!,:source]; etarget = al[!,:target];

julia> etype = map(x -> (startswith(x, "H") ? "admix" : "edge"), etarget)

julia> net1 = edgelist_to_net(esource, etarget, etype); writenewick(net1) # edge weights default to 0.5
"((t4:0.5,((t3:0.5,t1:0.5)i1:0.5,(t2:0.5)#H1:0.0::0.5)i2:0.5)i3:0.5,#H1:0.0::0.5)i4;"

julia> hardwiredclusterdistance(net0, net1, false) # false: as *un*rooted networks
0
```
"""
function edgelist_to_net(
    edgesource::AbstractVector,
    edgetarget::AbstractVector,
    edgetype::AbstractVector,
    edgeweight::Union{AbstractVector,Missing}=missing,
    weighttol = -1e-8
)
  nedges = length(edgesource)
  length(edgetarget) == nedges || error("edge source & target vectors should have same length")
  length(edgetype) == nedges || error("edge type vector should have length $nedges")
  if ismissing(edgeweight)
    edgeweight = repeat([0.5], nedges) # default edge lengths and γ = 0.5
  else
    length(edgeweight) == nedges || error("edge weight vector should have length $nedges")
    # check that edge weights are non-negative, with some tolerance
    for i in eachindex(edgeweight)
        edgeweight[i] < 0.0 || continue
        edgeweight[i] >= weighttol ||
            @error("negative edge weights $(edgeweight[i]) < tolerated $weighttol, for edge # i=$i")
        edgeweight[i] = 0.0
    end
  end
  nodenames = union(edgesource, edgetarget)
  # create all nodes
  nodevec = PhyloNetworks.Node[]
  nodename2index = Dict{eltype(nodenames),Int}()
  isrootnode = trues(length(nodenames)) # to find root later
  for (ni,nname) in enumerate(nodenames)
    nn = PhyloNetworks.Node(ni, true) # leaf by default
    nn.name = nname
    push!(nodevec, nn)
    push!(nodename2index, nname => ni)
  end
  # create all edges, and connects nodes & edges
  edgevec = PhyloNetworks.Edge[]
  ei = 0
  for (snn, tnn, et, ew) in zip(edgesource, edgetarget, edgetype, edgeweight)
    ei += 1
    sn = nodevec[nodename2index[snn]]
    tn = nodevec[nodename2index[tnn]]
    hyb = et == "admix"
    if hyb
        ee = PhyloNetworks.Edge(ei, 0.0)
        ee.hybrid = true
        ee.gamma = ew
        tn.hybrid = true
    else
        ee = PhyloNetworks.Edge(ei, ew)
    end
    sn.leaf = false # the source is not a leaf
    isrootnode[nodename2index[tnn]] = false # the target is not the root
    push!(ee.node, tn) # child first, because ischild1 is true by default
    push!(ee.node, sn)
    push!(edgevec, ee)
    push!(sn.edge, ee)
    push!(tn.edge, ee)
  end
  # rename hybrid nodes to Hi with i = i0,...,i0+h
  ihyb = 0
  for nn in nodevec # handle rare case when a tree node is already named "Hi"
    nn.hybrid && continue
    m = match(r"^H(\d+)$", nn.name)
    isnothing(m) && continue
    ihyb = max(ihyb, parse(Int, m.captures[1]))
  end
  for nn in nodevec
    nn.hybrid || continue
    ihyb += 1
    nn.name = "H$ihyb"
  end
  # for each hybrid node: normalize its parents' γs, find major parent
  for nn in nodevec
    nn.hybrid || continue
    hybparents = [ee for ee in edgevec if getchild(ee).name == nn.name]
    np = length(hybparents)
    np > 1 || error("hybrid node $nn with $np parents: $hybparents")
    gammasum = sum(ee.gamma for ee in hybparents)
    for ee in hybparents
        ee.gamma /= gammasum
        ee.ismajor = false
    end
    majorparentindex = argmax(ee.gamma for ee in hybparents)
    hybparents[majorparentindex].ismajor = true
  end
  # find the root
  sum(isrootnode) == 1 || error("found $(sum(isrootnode)) root(s)")
  rootindex = findfirst(isrootnode)
  # create the network
  net = PhyloNetworks.HybridNetwork(nodevec, edgevec)
  net.rooti = rootindex
  return net
end
