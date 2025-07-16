#=
Define networks of various complexities to fill the gap between g1 and g4.

g1_lg: deepcopy of fleg_pruned: large γ=0.48
g1_sg: same as g1_lg (h=1 reticulation) but with small γ=0.06

g2 with h=2 reticulations:
* one with γ=0.48 (into Denisovan if kept, into non-African otherwise)
* the other reticulation with γ=0.32.

- level 1: separate cycles, numbers give the cycle size for each reticulation
  * l1_c33, l1_c43, l1_c44, l1_c45
- level 2: numbers indicate the minimum cycle size for each reticulation,
  when the other reticulation is removed (by removing either hybrid edge)
  * l2_g33, l2_g46: galled and tree-child
  * l2_n04: with a hybrid ladder, so neither tree-child nor galled
  * l2_n44, l2_n46: tree-child but not galled

g3 with h=3 reticulations: both subnetworks of g4_lg with large γs,
g3_tc: tree-child (remove the parent hybrid in the stack)
g3_ntc: non tree-child (remove the ghost reticulation)
g3_l3_c54: (remove the ghost reticulation from g4_l4_c34 below)

g4_sg: deepcopy of fleg: one γ is small at 0.06
g4_lg: same as g4 (h=4 reticulations) but with large γs: 2 γ=0.48 and 2 γ=0.32
g4_l4_c34: from g4, moved origin of gene flow neanderthal_west → denisovan
=#

# load packages and define flegpopsizes & networks
#   with bottleneck: fleg (g4), fleg_pruned (g1)
#   and without: networkfleg_nobott, fleg_pruned_nobott, fleg_major 
include("fleg-net.jl")
const PN = PhyloNetworks

# edge numbers of minor hybrids:
[e.number => e.gamma for e in fleg.edge if !e.ismajor]
#=
  4 => 0.52  into Denisovan
  7 => 0.06  top of hybrid ladder into non-African
 11 => 0.32  bottom of hybrid ladder into non-African
 26 => 0.19  from ghost into base of modern humans
=#

"""
    findedge_bynumber(net, num)

edge in `net` whose number is `num`
"""
function findedge_bynumber(net::HybridNetwork, num::Integer)
    i = findfirst(x -> x.number == num, net.edge)
    isnothing(i) && return nothing # error("no edge numbered $num")
    return net.edge[i]
end

"""
    findnode_bynumber(net, num)

node in `net` whose number is `num`
"""
function findnode_bynumber(net::HybridNetwork, num::Integer)
    i = findfirst(x -> x.number == num, net.node)
    isnothing(i) && error("no node numbered $num")
    return net.node[i]
end

"""
    treeedgeγ1!(network)

Set all tree edges to have γ=1. Return the number of γs that needed to be fixed.
"""
function treeedgeγ1!(net)
    nfixed = 0
    for e in net.edge
        (e.hybrid || e.gamma == 1.0) && continue
        nfixed += 1
        e.gamma = 1.0
    end
  return nfixed
end

"""
    number_displayedtopologies(net)

number of distinct tree topologies displayed by `net`, ignoring degree-2 nodes.
Assumes that `net` has a single outgroup, to speed up the distance calculation
between trees assuming then as rooted (they all have the same outgroup if net
has a simple outgroup with no reticulation with the ingroup).
"""
function number_displayedtopologies(net)
  dtt = removedegree2nodes!.(displayedtrees(net, 0.0)) # suppresses the root too
  nd = length(dtt)
  for j in reverse(1:nd)
    for i in 1:(j-1)
      if hardwiredclusterdistance(dtt[j], dtt[i], true) == 0 # rooted=true okay
        deleteat!(dtt, j)
        break # don't try other i's!
      end
    end
  end
  return length(dtt)
end

"""
    blobsizes(net)

Non-trivial blob sizes. In a network with a single reticulation, this is the
number of edges in its unique cycle.
"""
function blobsizes(net)
  PN.process_biconnectedcomponents!(net)
  bs = [length(blob.edges) for blob in net.partition]
  filter!(!isequal(1), bs)
  return bs
end
"""
    subcyclesizes(net::HybridNetwork)
    print_subcyclesizes(net::HybridNetwork)

Cycle sizes in subnetworks: when keeping a single reticulation.
Example:

```julia
julia> subcyclesizes(fleg)
4-element Vector{Vector{Int64}}:
 [3, 4, 3, 3, 4, 4, 4, 4]
 [5, 4, 0, 0, 6, 5, 0, 0]
 [7, 7, 5, 6, 7, 6, 4, 4]
 [3, 4, 4, 5, 3, 4, 3, 4]

julia> print_subcyclesizes(fleg)
"34334444;54006500;77567644;34453434"

julia> map(v -> minimum(filter(x->x>0, v)), subcyclesizes(fleg))
4-element Vector{Int64}:
 3
 4
 4
 3

julia> sum(cs -> cs==3, map(v -> minimum(filter(x->x>0, v)), subcyclesizes(fleg)))
2
```
"""
function subcyclesizes(net::HybridNetwork)
  hybnode_num = [n.number for n in net.hybrid]
  nh = length(hybnode_num)
  hybedge_num = [ [getparentedge(n).number,getparentedgeminor(n).number]
    for n in net.hybrid]
  switchings = collect(Iterators.product([1:2 for _ in 1:(nh-1)]...))
  scs = Vector{Int}[]
  for i1 in 1:nh
    # keep hybedge_num[i1]. 2^{nh-1} ways to remove the others
    scs_i1 = Int[]
    he_num = deleteat!(copy(hybedge_num), i1)
    for I in switchings
      sunlet = deepcopy(net)
      for j2 in 1:(nh-1)
        en = he_num[j2][I[j2]]
        ee = findedge_bynumber(sunlet, en)
        if isnothing(ee)
          @info "i1=$i1, j2=$j2: edge $en was already removed (up a stack?)"
        else
          PN.deletehybridedge!(sunlet, ee)
        end
      end
      removedegree2nodes!(sunlet, false) # keeproot=false
      bs = blobsizes(sunlet)
      if length(bs) == 0 # hybedge_num[i1] was in a stack and removed
        push!(scs_i1, 0)
      elseif length(bs) == 1
        push!(scs_i1, bs[1])
      else
        @error "> 1 non-trivial blobs: bs=$bs. i1=$i1 will be ignored"
      end
    end
    push!(scs, scs_i1)
  end
  return scs
end

print_subcyclesizes(scs::Vector) = join(join.(scs, ""),";")
print_subcyclesizes(net::HybridNetwork) = print_subcyclesizes(subcyclesizes(net))

#--- networks with degree-2 nodes to keep using the map edge => population size

g1_lg = deepcopy(fleg_pruned);

g1_sg = deepcopy(fleg_pruned)
setgamma!(g1_sg.edge[findfirst(e -> e.gamma<0.5, g1_sg.edge)], 0.06)


g2_l1_c44 = deepcopy(fleg);
# below: delete hybrid edge with option nofuse=true to keep the degree-2 node
PN.deletehybridedge!(g2_l1_c44, findedge_bynumber(g2_l1_c44, 7),  true);
PN.deletehybridedge!(g2_l1_c44, findedge_bynumber(g2_l1_c44, 26), true);
[(number=e.number, γ=e.gamma) for e in g2_l1_c44. edge if e.gamma != 1 && !e.hybrid]
#  (number = 12, γ = 0.94)
#  (number = 22, γ = 0.81)
treeedgeγ1!(g2_l1_c44)  # 2

g2_l1_c45 = deepcopy(g2_l1_c44) 
αu = findedge_bynumber(g2_l1_c45, 16)
βu = getpartneredge(αu)
u = getchild(αu); u.number == 8
uv = findedge_bynumber(g2_l1_c45, 15)
v = getchild(uv); v.number == 11
vδ = findedge_bynumber(g2_l1_c45, 13)
savedlengths = [e.length for e in [αu,βu,uv,vδ]] # 387, 276, 568, 1976
nni!(αu, u, uv, v, vδ); # to detach uβ and graft it onto vδ
αu.length = savedlengths[1] + savedlengths[3] # 955
βu.length = savedlengths[2] + savedlengths[3] + 180 # 1024
uv.length = 180
vδ.length = savedlengths[4] - 180 # 1796
rotate!(g2_l1_c45, 11);
# cleanup namespace
αu = βu = u = uv = v = vδ = nothing
empty!(savedlengths)

g2_l1_c43 = deepcopy(g2_l1_c44);
PN.deletehybridedge!(g2_l1_c43, findedge_bynumber(g2_l1_c43, 11), true)
treeedgeγ1!(g2_l1_c43)  # 1
# code borrowed from addhybridedge! to add a hybrid edge
# from an existing node (-10) to an existing edge (19)
pn = findnode_bynumber(g2_l1_c43, -10)
ce = findedge_bynumber(g2_l1_c43,  19)
cn, pe = PN.breakedge!(ce, g2_l1_c43)
cn.hybrid = true
cn.name = "H4"
PN.pushHybrid!(g2_l1_c43, cn)
pe.hybrid = true
pe.gamma  = 0.68
flegpopsizes[pe.number] # new edge number: 28; so pop size Ne=27716
flegpopsizes[ce.number] # old edge 19: pop size Ne=33828, quite close: all good
he_length = sum(findedge_bynumber(g2_l1_c43, i).length for i in [22,21,28]) # 9329.5
he = PN.Edge(11, he_length, true, 0.32, false) # number, length, hybrid, γ, ismajor
PN.setNode!(he, [cn,pn])
PN.setEdge!(pn, he)
PN.setEdge!(cn, he)
PN.pushEdge!(g2_l1_c43, he)
PN.norootbelow!(ce)
for i in [-10,16] rotate!(g2_l1_c43, i); end
# root edge will be numbered 29. Was numbered 28 for fleg:
# define same population size for pop 29
push!(flegpopsizes, 29 => flegpopsizes[28])

g2_l1_c33 = deepcopy(g2_l1_c43);
ee = findedge_bynumber(g2_l1_c33, 23)
n18 = findnode_bynumber(g2_l1_c33, 18)
nm3 = findnode_bynumber(g2_l1_c33, -3)
PN.removeEdge!(n18, ee)
PN.removeNode!(n18, ee)
PN.setEdge!(nm3, ee)
PN.setNode!(ee, nm3)
directedges!(g2_l1_c33)
ee.length +=
  findedge_bynumber(g2_l1_c33, 24).length + # (-4 -> 18), 8960 generations
  findedge_bynumber(g2_l1_c33, 25).length   # (-3 -> -4), 5525 generations

g2_l2_g33 = deepcopy(fleg);
PN.deletehybridedge!(g2_l2_g33, findedge_bynumber(g2_l2_g33,7),  true);
PN.deletehybridedge!(g2_l2_g33, findedge_bynumber(g2_l2_g33,11), true);
setgamma!(findedge_bynumber(g2_l2_g33,26), 0.32) # fixit
[(number=e.number, γ=e.gamma) for e in g2_l2_g33.edge if e.gamma != 1 && !e.hybrid]
#  (number = 16, γ = 0.6799999999999999)
treeedgeγ1!(g2_l2_g33) # 1

g2_l2_g46 = deepcopy(fleg);
PN.deletehybridedge!(g2_l2_g46, findedge_bynumber(g2_l2_g46,12),  true);
PN.deletehybridedge!(g2_l2_g46, findedge_bynumber(g2_l2_g46,26), true);
[(number=e.number, γ=e.gamma) for e in g2_l2_g46.edge if e.gamma != 1 && !e.hybrid]
#   (number = 7, γ = 0.06)
#   (number = 22, γ = 0.81)
treeedgeγ1!(g2_l2_g46) # 2

g3 = deepcopy(fleg);
PN.deletehybridedge!(g3, findedge_bynumber(g3,4),  true);
treeedgeγ1!(g3) # 1

g2_l2_n04 = deepcopy(g3);
PN.deletehybridedge!(g2_l2_n04, findedge_bynumber(g2_l2_n04,26),  true);
setgamma!(findedge_bynumber(g2_l2_n04,11), 0.48)
setgamma!(findedge_bynumber(g2_l2_n04,7),  0.32)
treeedgeγ1!(g2_l2_n04) # 1

g2_l2_n44 = deepcopy(g3);
PN.deletehybridedge!(g2_l2_n44, findedge_bynumber(g2_l2_n44,7),  true);
setgamma!(findedge_bynumber(g2_l2_n44,11), 0.48)
setgamma!(findedge_bynumber(g2_l2_n44,26), 0.32)
treeedgeγ1!(g2_l2_n44) # 1

g2_l2_n46 = deepcopy(g3);
PN.deletehybridedge!(g2_l2_n46, findedge_bynumber(g2_l2_n46,12),  true);
setgamma!(findedge_bynumber(g2_l2_n46,11), 0.48)
setgamma!(findedge_bynumber(g2_l2_n46,26), 0.32)
treeedgeγ1!(g2_l2_n46) # 1

g3=nothing # was introduced to define g2_l2_n*, no longer needed

g1_l1_c3 = deepcopy(g2_l1_c43);
PN.deletehybridedge!(g1_l1_c3, findedge_bynumber(g1_l1_c3,4),  true);
treeedgeγ1!(g1_l1_c3) # 1

g4_sg = deepcopy(fleg);

g4_lg = deepcopy(fleg);
setgamma!(findedge_bynumber(g4_lg,11), 0.48)
setgamma!(findedge_bynumber(g4_lg, 7), 0.32)
setgamma!(findedge_bynumber(g4_lg,26), 0.32)

g3_tc = deepcopy(g4_lg);
PN.deletehybridedge!(g3_tc, findedge_bynumber(g3_tc,7),  true);
treeedgeγ1!(g3_tc) # 1

g3_ntc = deepcopy(g4_lg);
PN.deletehybridedge!(g3_ntc, findedge_bynumber(g3_ntc,22),  true);
treeedgeγ1!(g3_ntc) # 1

g4_l4_c34 = deepcopy(g4_lg);
αu = findedge_bynumber(g4_l4_c34, 5)
βu = findedge_bynumber(g4_l4_c34, 8)
u = getparent(αu); u.number == 7
uv = findedge_bynumber(g4_l4_c34, 9)
v = getparent(uv); v.number == -7
vδ = findedge_bynumber(g4_l4_c34, 10)
vγ = findedge_bynumber(g4_l4_c34, 4)
savedlengths = [e.length for e in [αu,βu,uv,vδ,vγ]] # 223, 1193, 11933, 4309, 8635
nni!(αu, u, uv, v, vδ); # to detach uβ and graft it onto vδ
vδ.length = savedlengths[4] + 3000 # 7309
βu.length = savedlengths[2] + (savedlengths[3] - 3000) # 10126
uv.length = savedlengths[5] - 3000 # 5635
vγ.length = 0
αu.length = savedlengths[1] + savedlengths[3] - savedlengths[5] # 3521
rotate!(g4_l4_c34, 7);
# plot(g4_l4_c34, useedgelength=true, showedgenumber=true, xlim=[220_000, 242_000]);
istimeconsistent(g4_l4_c34) # should be true
# cleanup namespace
αu = βu = u = uv = v = vδ = nothing
empty!(savedlengths)

#= move recipient of gene flow deep ghost → africa_west
#  but: fails to get all minimum cycles of 4+ nodes
g4_l4_c43 = deepcopy(g4_l4_c34);
αu = findedge_bynumber(g4_l4_c43, 22)
βu = getpartneredge(αu)
u = getchild(αu); u.number == 17
uv = findedge_bynumber(g4_l4_c43, 21)
v = getchild(uv); v.number == 16
vδ = findedge_bynumber(g4_l4_c43, 20)
vγ = findedge_bynumber(g4_l4_c43, 19)
savedlengths = [e.length for e in [αu,βu,uv,vδ,vγ]] # 2523, 23079, 5347, 6314, 2919
nni!(αu, u, uv, v, vδ); # to detach uβ and graft it onto vδ
αu.length = savedlengths[3]
βu.length = savedlengths[2] + savedlengths[3]
uv.length = savedlengths[1]
vδ.length = savedlengths[4]
vγ.length = savedlengths[5] + savedlengths[1]
istimeconsistent(g4_l4_c43) # should be true
# cleanup namespace
αu = βu = u = uv = v = vδ = nothing
empty!(savedlengths)
# bummer: new 3-cycles appeared, at a hybrid that used to have 4+ cycles:
print_subcyclesizes(g4_l4_c43) # 45445555;44005500;67566644;34453434

# possibly, to get g4_l4_c44: add operation done to g3_ntc
=#

g3_l3_c54 = deepcopy(g4_l4_c34);
PN.deletehybridedge!(g3_l3_c54, findedge_bynumber(g3_l3_c54,26),  true);
treeedgeγ1!(g3_l3_c54) # 1

#------- define networks without bottleneck ------#

population_net_names = [
  :g1_lg, :g1_sg, :g1_l1_c3,
  :g2_l1_c33, :g2_l1_c43, :g2_l1_c44, :g2_l1_c45,
  :g2_l2_g33, :g2_l2_g46, :g2_l2_n04, :g2_l2_n44, :g2_l2_n46,
  :g3_tc, :g3_ntc, :g3_l3_c54,
  :g4_l4_c34, :g4_lg, :g4_sg
]
for netname in population_net_names
    net = eval(netname)
    global tmpnet
    tmpnet = deepcopy(net)
    removedegree2nodes!(tmpnet, true)
    tmpname = Symbol("$netname" * "_nobot")
    eval(quote
        $tmpname = tmpnet
    end)
end

#=
old2newnames = Dict( # to more easily describe information contained in names
  "g2-l2-n04" => "g2-l2-n54", # 'n' for neither galled nor tree-child
  "g2-l2-n44" => "g2-l2-t44", # 't' for tree-child (only)
  "g2-l2-n46" => "g2-l2-t46",
  "g3-tc"  => "g3-t334",
  "g3-ntc" => "g3-n345",
  "g3-l3-c54" => "g3-l3-c545",
  "g4-l4-c34" => "g4-l4-c3444",
  "g1-sg" => "g1-lowγ",
  "g4-lg" => "g4-highγ",
);
R"cairo_pdf"("figures/SMfig_popnets.pdf", height=11, width=11);
R"layout"([1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]);
R"par"(mar=[0,0,0,0]);
for (ii, netname) in enumerate(population_net_names)
  netname in [:g1_lg,:g4_sg] && continue # skip: topology repeat
  net = eval(netname)
  for n in net.leaf # shorten taxon names
    n.name = replace(n.name, "chimp" => "C", "denisovan" => "D", "neanderthal" => "N",
      "non_africa" => "nA", "africa" => "A", "_east" => "e", "_west" => "w")
  end
  plot(net, showgamma=true, showedgelength=false, tipoffset=0.1);
  netlabel = replace(string(netname), "_"=>"-")
  for (old,new) in old2newnames netlabel = replace(netlabel, old => new); end
  R"mtext"(text=netlabel, line=-1.3, cex=0.8, adj=0.1)
  R"mtext"(text="x=$(nt_popnetinfo.num_mincyclesize3[ii])", line=-2.3, cex=0.8, adj=.09)
end
R"dev.off"();

# simple plot: no names, no leaf labels
R"cairo_pdf"("../figures/fig_popnets_simple.pdf", height=3, width=4);
R"layout"([1 2 3 4; 5 6 7 8; 9 10 11 12]);
R"par"(mar=[0,0,0,0]);
for netname in population_net_names
  netname == :g1_l1_c3 && continue # skip plotting this one for now
  net = eval(netname)
  plot(net, showtiplabel=false, majorhybridedgecolor="black", arrowlen=0.03);
end
R"dev.off"();

# g2 c43 and g33 only: to highlight level and minimum cycle size
nn_c43 = DataFrame(nodenumber=[19,16,-10, 3,-7,18,-4],
  col=["orangered","orangered","orangered", "blue","blue","blue","blue"],
  lab=repeat([""],inner=7))
nn_g33 = DataFrame(nodenumber=[17,-3, -4,18, 3,-7],
  col=["orangered","orangered","tan","tan","blue","blue"],
  lab=repeat([""],inner=6))
R"cairo_pdf"("../figures/fig_popnets_c43_g33.pdf", height=3, width=6);
R"layout"([1 2]);
R"par"(mar=[0,0,0,0]);
for netname in [:g2_l1_c43_nobot, :g2_l2_g33_nobot]
  net = eval(netname)
  nn = (netname==:g2_l1_c43_nobot ? nn_c43 : nn_g33)
  p = plot(net, xlim=(0.8,9), majorhybridedgecolor="deepskyblue",
    showtiplabel=false, nodelabel=nn[:,[:nodenumber, :lab]]);
  ndf = p[:node_data]
  o = indexin(nn.nodenumber, parse.(Int,ndf.num))
  R"points"(x=ndf.x[o], y=ndf.y[o], pch=16, col=nn.col, cex=1.5);
end
R"dev.off"();

# figure to explain minimum cycle size, :g2_l2_g33 to illustrate
using DataFrames
nn_g33 = DataFrame(nodenumber=[17,-3, 18,-4, 3,-7],
  col=["orangered4","orangered4","orangered","dodgerblue","dodgerblue4","dodgerblue4"],
  lab=repeat([""],inner=6))
R"cairo_pdf"("figures/fig_mincyclesize_g33.pdf", height=2, width=5.5);
R"layout"([1 2 3]);
R"par"(mar=[0,0,0,0]);
ecd = Dict(3=>"deepskyblue", 4=>"deepskyblue", 22=>"orangered", 26=>"orangered")
pargs = (xlim=(0.8,9), showtiplabel=false, edgecolor=ecd,)
p = plot(g2_l2_g33_nobot; nodelabel=nn_g33[:,[:nodenumber, :lab]], pargs...);
ndf = p[:node_data]
o = indexin(nn_g33.nodenumber, parse.(Int,ndf.num))
ndf = ndf[o,:]
R"points"(x=ndf.x, y=ndf.y, pch=16, col=nn_g33.col, cex=1.5);
ecd[26] = "white"
plot(g2_l2_g33_nobot; preorder=false, pargs...);
R"points"(x=ndf.x[3:6], y=ndf.y[3:6], pch=16, col=nn_g33.col[3:6], cex=1.5);
ecd[26] = "orangered"; ecd[22] = "white"
plot(g2_l2_g33_nobot; preorder=false, pargs...);
R"points"(x=ndf.x[4:6], y=ndf.y[4:6], pch=16, col=nn_g33.col[4:6], cex=1.5);
R"dev.off"();
# then edited with illustrator to remove extra segments

=#

nt_popnetinfo = ( # NamedTuple. can convert to a DataFrame later
  net=String[],
  trueh=Int[],
  level=Int[],
  num_dtrees=Int[],
  gamma=String[],
  num_mincyclesize3=Int[],
  cyclesizes=String[],
)
for netname in population_net_names
  net = eval(netname)
  netname_dash = replace(string(netname), "_"=>"-")
  smallγ = any(e.gamma < 0.10 for e in net.edge)
  h = parse(Int, match(r"g(\d)-", netname_dash).captures[1])
  level = h # correct if h=1, 3 or 4; not necessarily if h=2
  if level==2
    level = parse(Int, match(r"g2-l(\d)-", netname_dash).captures[1])
  end
  ndt = number_displayedtopologies(net)
  scs = subcyclesizes(net)
  # number of hybridizations whose minimum cycle size is 3
  #   after filtering out 0-cycles: when the reticulation is non-visible
  n3c = sum(cs -> cs==3, map(v -> minimum(filter(x->x>0, v)), scs))
  push!(nt_popnetinfo.net,   netname_dash)
  push!(nt_popnetinfo.level, level)
  push!(nt_popnetinfo.num_dtrees, ndt)
  push!(nt_popnetinfo.num_mincyclesize3, n3c)
  push!(nt_popnetinfo.cyclesizes, print_subcyclesizes(scs))
  push!(nt_popnetinfo.gamma, (smallγ ? "small" : "large"))
  push!(nt_popnetinfo.trueh, h)
end
#=
julia> using CSV, DataFrames

julia> DataFrame(nt_popnetinfo)
18×7 DataFrame
 Row │ net        trueh  level  num_dtrees  gamma   num_mincyclesize3  cyclesizes                        
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ g1-lg          1      1           2  large                   0  4
   2 │ g1-sg          1      1           2  small                   0  4
   3 │ g1-l1-c3       1      1           1  large                   1  3
   4 │ g2-l1-c33      2      1           1  large                   2  33;33
   5 │ g2-l1-c43      2      1           2  large                   1  44;33
   6 │ g2-l1-c44      2      1           4  large                   0  44;44
   7 │ g2-l1-c45      2      1           4  large                   0  44;55
   8 │ g2-l2-g33      2      2           2  large                   2  34;34
   9 │ g2-l2-g46      2      2           4  large                   0  44;76
  10 │ g2-l2-n04      2      2           3  large                   0  50;64
  11 │ g2-l2-n44      2      2           4  large                   0  64;54
  12 │ g2-l2-n46      2      2           4  large                   0  76;44
  13 │ g3-tc          3      3           6  large                   2  4344;5644;4534
  14 │ g3-ntc         3      3           4  large                   1  3433;5400;7756
  15 │ g3-l3-c54      3      3           6  large                   0  5555;5500;6644
  16 │ g4-l4-c34      4      4          10  large                   1  45445555;44005500;67566644;34453434
  17 │ g4-lg          4      4           8  large                   2  34334444;54006500;77567644;34453434
  18 │ g4-sg          4      4           8  small                   2  34334444;54006500;77567644;34453434

julia> CSV.write("pop-net-properties.csv", nt_popnetinfo)
=#

#=
export these graphs into the "edgelist" format that can be used by admixtools
then save them in csv file that can be read without re-running this julia code
also:
- save their edge lengths in coalescent units
- remove bottleneck
- but add degree-2 node at the end of any hybrid edge of positive length
done with code below, mostly commented out because output csv files are tracked
=#

#=
using DataFrames, CSV
# load definition of getedgelist_coalunits!:
include("../scripts/interop_admixtools.jl")
for netname in population_net_names
  net = eval(netname)
  netname_dash = replace(string(netname), "_"=>"-")
  filename = joinpath("input/pop-net_edgelist", netname_dash * "_cu.csv")
  getedgelist_coalunits!(deepcopy(net), popsizedict=flegpopsizes, csvfile=filename)
end
=#

@info """now defined:
functions: findedge_bynumber, findnode_bynumber, treeedgeγ1!,
number_displayedtopologies, blobsizes, subcyclesizes, print_subcyclesizes

$(length(population_net_names)) population networks, all including the bottleneck node:
- g1_lg = fleg_pruned and g1_sg with small γ
  g1_l1_c3: 3-cycle with semi-large γ (0.32) as in g2_l1_c43
- h=2, with semi-large γs in (0.32 and 0.48):
  g2_l1_c33, g2_l1_c43, g2_l1_c44, g2_l1_c45 (level-1)
  g2_l2_g33, g2_l2_g46 (level-2, galled, tree-child)
  g2_l2_n04, g2_l2_n44, g2_l2_n46 (level-2, not tree-child or not galled)
- h=3, with semi-large γs in (0.32 and 0.48):
  g3_tc (tree-child), g3_ntc (not tree-child), g3-l3-c54 (all cycles 4+ nodes)
- g4_sg = fleg, g4_lg with semi-large γs, g4-l4-c34

networks with all degree-2 nodes removed: g1_sg_nobot, g2_l1_c44_nobot etc.
population_net_names : vector of symbols, containing the network's names
"""
