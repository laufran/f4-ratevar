using PhyloPlots
using PhyloNetworks
using RCall

#=
source of network & edge parameters (pop. size and lengths in # of generations):
Flegontov et al. (2023) preprint: https://doi.org/10.1101/2023.01.22.525077
from their Fig. 2a, and Supp. Table 13 for population sizes not shown in Fig. 2a.

taxon names ending with _1 (or _2) were changed to end with _west (or _east)
to ease the separation between population name and individual number.

final paper in PLOS Genetics 19(9):e1010931
Flegontov et al. (2023) https://doi.org/10.1371/journal.pgen.1010931
from which this network is Fig. 3a, and other parameters in S13 Table.
"The former admixture graph (Fig 3A) reproduces some known features of the genetic
history of anatomically modern and archaic humans, but differs in other respects
from the widely accepted model [32,53]"
[32] Prufer et al. 2014 https://doi.org/10.1038/nature12886
[53] Durvasula & Sankararaman 2020 https://doi.org/10.1126/sciadv.aax5097
=#

fleg_str = """
(chimp:240000,(
  
((denisovan:5611)#H2:21904,((#H2:8635::0.52[&neanderthal2denisovan],(neanderthal_west:223,(neanderthal_east:1120,#H3:0::0.06):1193)Ne0:11933):4309,
   
((#H4:276::0.32[&neanderthals2nonafrican])#H3:11364,((((((non_africa_west:1976,non_africa_east:1976)na0:568)#H4:387)bottleneck:464,africa_east:3395)a0:2919,africa_west:6314)AMH:5347)#H1:2523):6071)NeAMH:8960
  ):5525,
  #H1:23079::0.19[&ghostarchaic]
):205260);
"""
fleg = readnewick(fleg_str)
fleg.edge[3].ismajor = true # to restore viz in fig2a of Flegontov
fleg.edge[4].ismajor = false # despite γ = 52% for minor edge
#=
R"pdf"("flegontov2023-fig2a.pdf", height=5, width=8)
R"par"(mar=[0,0,0,0]);
plot(fleg, shownodelabel=false, xlim=[0,14], showgamma=true,
     edgelabel = DataFrame(
        num = [e.number for e in fleg.edge],
        lab = [string(Int(e.length)) for e in fleg.edge],
     ), edgecex=0.5, tipoffset=0.2);
R"dev.off"();
=#

#=
# visual check of edge lengths -- but too rough
hominin = deepcopy(fleg); deleteleaf!(hominin, "chimp");
plot(hominin, shownodelabel=true, useedgelength=true);
=#

# population sizes in number of diploid individuals, from fig. 2a and Supp. Table 13
flegpopsizes = Dict{Int64, Int64}(
    1 => 1000, #chimpanzee (outgroup)
    2 => 16758, #Denisovan (d)
    3 => 16758, #hybrid edge w/ gamma = 48%
    4 => 7311, #hybrid edge w/ gamma=52%
    5 => 14399, #neanderthal_west
    6 => 14399, #neanderthal_east
    7 => 2820, #hybrid edge w/ gamma = 6%
    8 => 14399,
    9 => 8145, #ancestral Neanderthal population
    10 => 8145, 
    11 => 2544, #hybrid edge w/ gamma=32%
    12 => 2544, #hybrid edge w/ gamma=94% 
    13 => 35744, #non_africa_west
    14 => 14763, #non_africa_east
    15 => 8821, #ancestral non-African population (na_w + na_e)
    16 => 8821, #major edge (68%) to 11
    17 => 1506, #out-of-Africa bottleneck
    18 => 46139, #africa_east
    19 => 16914, #a0 (ancestor of a_e + na_w + na_e)
    20 => 44541, #africa_west
    21 => 222379, #AMH (a_w + a_e + na_w + na_e)
    22 => 222379, #hybrid edge w/ gamma=81%
    23 => 222379,
    24 => 86161, #Neanderthal + AMH ancestral population
    25 => 35414,
    26 => 11661, #hybrid edge w/ gamma=19%
    27 => 35414, #super-archaic population
    28 => 13858 #root population (chimpanzee + archaic + AMH)
)

# convert from diploid to haploid population sizes
for k in keys(flegpopsizes)
  flegpopsizes[k] *= 2
end

# remove degree 2 nodes from fleg unpruned (remove bottleneck)
# -> need this for tree/net distance calculations
fleg_nobott = deepcopy(fleg)
removedegree2nodes!(fleg_nobott, true)

#get pruned version of tree to use as starting tree
fleg_pruned = deepcopy(fleg)
deletehybridthreshold!(fleg_pruned, 0.35, true)
# nofuse = true: to have same edges as in fleg and re-use flegpopsizes

# remove degree 2 nodes- otherwise keeps nodes from pruned hybrid events, despite removing edges
fleg_pruned_nobott = deepcopy(fleg_pruned)
removedegree2nodes!(fleg_pruned_nobott, true)
# PhyloPlots.plot(fleg_pruned, tipoffset=0.1, showedgenumber=true, shownodenumber=true, showgamma=true);

# major tree from Flegontov network
fleg_major = majortree(fleg)
removedegree2nodes!(fleg_major, true)
#plot(fleg_major, showedgenumber = true) #confirm bottleneck is removed

#=
# plot network with edge lengths in generations, and population sizes labeled
plot(fleg, tipoffset=0.1, showedgelength=true, edgelabelcolor="red4",
          edgelabel=DataFrame(n=[e.number for e in fleg.edge],
                              l=[flegpopsizes[e.number] for e in fleg.edge]));
R"text"(x=1, y=2.5, flegpopsizes[28], adj=1, col="red4");
R"mtext"("red: Ne values", side=1, line=-1.5, col="red4");
R"mtext"("black: # gens", side=1, line=-0.5);

plot(fleg_pruned, tipoffset=0.1, showedgelength=true, edgelabelcolor="red4",
          edgelabel=DataFrame(n=[e.number for e in fleg_pruned.edge],
                              l=[flegpopsizes[e.number] for e in fleg_pruned.edge]));
R"text"(x=1, y=2.5, flegpopsizes[28], adj=1, col="red4");
R"mtext"("red: Ne values", side=1, line=-1.5, col="red4");
R"mtext"("black: # gens", side=1, line=-0.5);
=#

@info """now defined:
- flegpopsizes
- fleg (h=4), fleg_pruned (h=1)
- fleg_nobott (h=4), fleg_pruned_nobott (h=1), fleg_major (h=0): no bottleneck
"""
