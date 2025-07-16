#=
script to fix a bug in the calculation of the number of displayed
near major trees between each top graph and the true graph.

previously, we did not update the `num_nearmajortrees` if the 
`num_nearmajortrees < 2` initially, before changing the gamma threshold.
(addition of line 62).

Workflow:
- read CSV file with all topgraphs: adjust inputcsv file name below
- recalculate values in column ndisplayed
- save full data frame with updated column in CSV file, ending in "-june11.csv"
=#

# adjust below, to use the most recent after update
inputcsv = "output/concat_g2topgraphs_filtered-june15.csv"
outputcsv = replace(inputcsv, r"(-\w+)?\.csv" => "-june15_ngraphsfix.csv")

using PhyloNetworks
using CSV
using DataFrames
include("../input/pop-net.jl")

function ndisplayed_nearmajortrees(net, trees2find=nearmajortree, nmt=num_nearmajortrees)
  dtrees = removedegree2nodes!.(displayedtrees(net, 0.0))
  notfound = collect(eachindex(trees2find))
  for tree in dtrees
    isempty(notfound) && break
    for (ii,i) in enumerate(notfound)
      if hardwiredclusterdistance(tree, trees2find[i], false) == 0
        deleteat!(notfound, ii) # no need to look for trees2find[i] in later displayed trees
        break # 'tree' cannot match other near-major trees
      end
    end
  end
  return nmt - length(notfound)
end

#load in CSV
topgraphs = DataFrame(CSV.File(inputcsv))

#extract current ndisplayed- some with bug, some without
ndisplayed_old = topgraphs[:, :ndisplayed]
ndisplayed_new = deepcopy(ndisplayed_old)

#extract net names simulated under
nets = topgraphs[:, :net]

replace!(nets, "g2-l1-c4" =>"g2-l1-c44",
               "g2-l1-c5" =>"g2-l1-c45",
               "g2-l2-c4" =>"g2-l2-g33",
               "g2-l2-c6" =>"g2-l2-g46",
               "g2-l2-ntc"=>"g2-l2-n04",
               "g2-l2-ng4"=>"g2-l2-n44",
               "g2-l2-ng6"=>"g2-l2-n46")

#extract the newick cols to iterate over (not row by row in CSV)
newicks = topgraphs[:, :newick]

for idx in eachindex(nets)
    truenetname = nets[idx]
    newickstr = newicks[idx]

    netname = replace(truenetname, "-" => "_")
    netname = Symbol(netname * "_nobot")
    truenet = eval(netname)

    nearmajortree = displayedtrees(truenet, 0.45)
    num_nearmajortrees = length(nearmajortree)

    if num_nearmajortrees < 2
        if truenet.numhybrids == 1 # then use the 2 displayed tree, even if low γ
            nearmajortree = displayedtrees(truenet, 0.0)
            num_nearmajortrees = length(nearmajortree)
            newick = readTopology(newickstr)
            ndisplayed_new[idx] = ndisplayed_nearmajortrees(newick, nearmajortree, num_nearmajortrees)
        else
            @error "$num_nearmajortrees near-major trees in truenet, instead of 2"
        end
    end
end

#delete the old hw tree dist col
select!(topgraphs, Not([:ndisplayed]));
#replace with the new col
insertcols!(topgraphs, 10, :ndisplayed => ndisplayed_new)

#rewrite CSV
CSV.write(outputcsv, topgraphs)