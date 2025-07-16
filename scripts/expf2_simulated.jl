#= averaged observed f2s:
- as calculated by admixtools
- on the simulated data:
  * under our simulation model
  * and filtering to only keep biallelic sites

use `avgf2.csv`, 1 per replicate; average them across 100 replicates from each
scenario (combination of parameters). Warning: populations can be listed in
different orders across different replicates.

example use from output folder (runs in ~ 30 sec.):

```julia
include("../scripts/expf2_simulated.jl")
df = avgf2([r"fleg"])
CSV.write("averageobservedf2.csv", df);
```

=#

# load script that defines parse_folder2parameters:
include("probability_012mutations.jl")

pops = sort!(
    ["chimp", "denisovan", "neanderthal_west", "neanderthal_east",
    "non_africa_west", "non_africa_east", "africa_east", "africa_west"])
# before sorting: same order as in tiplabels(fleg)
poppair_index = [i=>j for i in 1:length(pops) for j in (i+1):length(pops)]
poppair_names = [pops[p.first] => pops[p.second] for p in poppair_index]

"""
    avgf2_onerep(folder; asmatrix::Bool=true)

f2 values from one file `avgf2.csv` (averaged across all sites in the alignment)
in `folder`, between populations `pops`.

By default, f2 values are returned in a vector, each value corresponding to one
pair. In this vector, pairs of populations are listed as in `poppair_names`,
whose indices in `pops` are assumed to be `poppair_index`.
These variables are external to the function.

With option `asmatrix=false`, a matrix is returned, with populations listed
in rows and columns in the same order as in `pops`.

# example

```julia
julia> tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate/001";

julia> avgf2_onerep(tmp; asmatrix=true)
8×8 Matrix{Float64}:
 0.0     0.2237  0.6403  0.2334  0.2295  0.2286  0.2062  0.2013
 0.2237  0.0     0.6318  0.2269  0.2242  0.2241  0.2237  0.215
 0.6403  0.6318  0.0     0.6289  0.6366  0.6345  0.6339  0.6336
 0.2334  0.2269  0.6289  0.0     0.2045  0.2034  0.2322  0.2249
 0.2295  0.2242  0.6366  0.2045  0.0     0.0913  0.2263  0.224
 0.2286  0.2241  0.6345  0.2034  0.0913  0.0     0.2282  0.2265
 0.2062  0.2237  0.6339  0.2322  0.2263  0.2282  0.0     0.1665
 0.2013  0.215   0.6336  0.2249  0.224   0.2265  0.1665  0.0

julia> avgf2_onerep(tmp) |> print
[0.2237, 0.6403, 0.2334, 0.2295, 0.2286, 0.2062, 0.2013, 0.6318, 0.2269, 0.2242, 0.2241, 0.2237, 0.215, 0.6289, 0.6366, 0.6345, 0.6339, 0.6336, 0.2045, 0.2034, 0.2322, 0.2249, 0.0913, 0.2263, 0.224, 0.2282, 0.2265, 0.1665]

```
"""
function avgf2_onerep(dir; asmatrix::Bool=false)
    file = joinpath(dir, "avgf2.csv")
    df = CSV.read(file, DataFrame)
    names(df)[2:end] == df.Column1 || # sanity check
        error("populations ordered differently in rows and columns. dir = $dir")
    select!(df, Not(:Column1))
    # Set(names(df)) == Set(pops) || error("unexpected set of populations. dir = $dir")
    o = indexin(pops, names(df))
    names(df)[o] == pops || error("unexpected set of populations. dir = $dir")
    np = length(pops)
    if asmatrix
        res = zeros(np, np)
        for ((i,j)) in poppair_index
            res[i,j] = res[j,i] = df[o[i],o[j]]
        end
    else
        res = zeros(np*(np-1)÷2)
        for (ij,(i,j)) in enumerate(poppair_index)
            res[ij] = df[o[i],o[j]]
        end
    end
    return res
end

"""
    avgf2_matrix(folder)

Matrix of averaged observed f2 values, averaged across all 100 replicates in the
`folder` from one scenario (1 combination of parameter values).

# example
```julia
julia> nind = [1,2,10];

julia> inname = ["fleg_$i-indls_118000-genes_100000-biall_0.0-lindist_1.25e-8-subrate" for i in nind];

julia> f2s = [avgf2_matrix(dir) for dir in inname];

julia> outname = joinpath.(inname, "avgf2.csv");

julia> for i in 1:3 CSV.write(outname[i], DataFrame(f2s[i], pops)); end
```
"""
function avgf2_matrix(dir)
    repdir = filter(contains(r"\d+"), readdir(dir))
    f2s = Matrix{Float64}[]
    for repstr in repdir
        # irep = parse(Int, repstr)
        repdir = joinpath(dir, repstr)
        push!(f2s, avgf2_onerep(repdir; asmatrix=true))
    end
    nreps = length(f2s)
    nreps == 100 || @warn "weird: $nreps replicates for $dir"
    f2mean = round.(mean(f2s), digits=10)
    return f2mean
end

"""
    avgf2(paramfilters::Vector{Regex}=Regex[])

DataFrame of averaged observed f2 values with one row per scenario (combination
of parameter values) that pass the filter(s). where f2s are averaged across the
100 replicates from each scenario.

# example

```julia
julia> paramf = ["fleg-pruned_", "_1.25e-7-subrate", "100000-biall"];

julia> df = avgf2(Regex.(paramf))
[ Info: starting: 15 parameter combination folders
[ Info: starting: fleg-pruned_1-indls_118000-genes_100000-biall_0.0-lindist_1.25e-7-subrate
...
15×34 DataFrame
 Row │ net          nind   nbiallelic  lindist  mutrate  nreps  africa_east__africa_west  africa_east__chimp  africa_east__denisovan  ⋯
     │ String       Int64  Int64       Float64  Float64  Int64  Float64                   Float64             Float64                 ⋯
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ fleg-pruned      1      100000     0.0   1.25e-7    100                 0.221667             0.637596               0.229965   ⋯
   2 │ fleg-pruned      1      100000     0.15  1.25e-7    100                 0.220999             0.638926               0.231622
   3 │ fleg-pruned      1      100000     0.3   1.25e-7    100                 0.223164             0.633722               0.229787
   4 │ fleg-pruned      1      100000     0.5   1.25e-7    100                 0.227939             0.621016               0.239792
   5 │ fleg-pruned      1      100000     0.7   1.25e-7    100                 0.21677              0.630965               0.227873   ⋯
   6 │ fleg-pruned     10      100000     0.0   1.25e-7    100                 0.0238216            0.319044               0.0436516
   7 │ fleg-pruned     10      100000     0.15  1.25e-7    100                 0.0238238            0.318101               0.0434364
   8 │ fleg-pruned     10      100000     0.3   1.25e-7    100                 0.0236335            0.309766               0.0435832
   9 │ fleg-pruned     10      100000     0.5   1.25e-7    100                 0.0234043            0.310759               0.0428828  ⋯
  10 │ fleg-pruned     10      100000     0.7   1.25e-7    100                 0.0236782            0.306899               0.0432462
  11 │ fleg-pruned      2      100000     0.0   1.25e-7    100                 0.10393              0.500253               0.123126
  12 │ fleg-pruned      2      100000     0.15  1.25e-7    100                 0.104179             0.500703               0.122914
  13 │ fleg-pruned      2      100000     0.3   1.25e-7    100                 0.104759             0.491447               0.123891   ⋯
  14 │ fleg-pruned      2      100000     0.5   1.25e-7    100                 0.108012             0.481772               0.126475
  15 │ fleg-pruned      2      100000     0.7   1.25e-7    100                 0.102969             0.486341               0.11953
                                                                                                                     25 columns omitted

julia> filename = "averageobservedf2_" * join(replace.(paramf, "_" => ""), '_') * ".csv";
```
"""
function avgf2(
    paramfilters::Vector{Regex}=Regex[],
)
    df = DataFrame(net=String[],
            nind=Int[], nbiallelic=Int[], lindist=Float64[], mutrate=Float64[],
            nreps=Int[])
    for ((pop1, pop2)) in poppair_names
        colname = Symbol(pop1 * "__" * pop2)
        df[!,colname] = Float64[]
    end
    paramset_dirs = filter(contains(r"\w+_\d.+"), readdir())
    for f in paramfilters
        filter!(contains(f), paramset_dirs)
    end
    @info "starting: $(length(paramset_dirs)) parameter combination folders"
    for dir in paramset_dirs
        avgf2!(df, dir)
    end
    return df
end
function avgf2!(df::DataFrame, dir::AbstractString)
    net, nind, nbiallelic, lindist, mutrate = parse_folder2parameters(dir)
    repdir = filter(contains(r"\d+"), readdir(dir))
    @info "starting: $dir"
    # starttime = time()
    f2s = Vector{Float64}[]
    for repstr in repdir
        # irep = parse(Int, repstr)
        repdir = joinpath(dir, repstr)
        push!(f2s, avgf2_onerep(repdir))
    end
    nreps = length(f2s)
    nreps == 100 || @warn "weird: $nreps replicates for $dir"
    f2mean = round.(mean(f2s), digits=10)
    push!(df, (net,nind,nbiallelic,lindist,mutrate, nreps, f2mean...))
    # endtime = time()
    # elapsed = round((endtime - starttime) / 60, digits=3)
    # @info "$elapsed mins for $dir"
    GC.gc()
    return df
end

#= figure: show networks that fit the avg f2 (and expected f2) data
   as perfectly as g4.

These are: g4-l4-c34 and graphs in file
`output/g14_topgraphs_byrep_scoretol0.05_correctToB_bestscenario.csv`
on rows 20,28,19,2,30 (after the header).
See expf2_findgraphs.{qmd,r} for how these networks were identified to fit
near perfectly (up to numerical precision) and better than g4.

The code below was run interactively, but commented out in `if false` block
to skip when including the file.
=#

if false # to skip when doing include("scripts/expf2_simulated.jl")

include("input/pop-net.jl")
using CSV, DataFrames
csvfile = "output/g14_topgraphs_byrep_scoretol0.05_correctToB_bestscenario.csv"
df = CSV.read(csvfile, DataFrame)
rotate!(fleg_nobott, 14);
rotate!(g4_l4_c34_nobot, 14);
netrows = [2,19,28,20,30]
altg = [readnewick(df.newick[i]) for i in netrows]
altr = [df.irep[i] for i in netrows]
# shorten names
for net in (fleg_nobott, g4_l4_c34_nobot, altg...), n in net.leaf
  n.name = replace(n.name,
    "non_africa" => "nA", "africa" => "A", "_east" => "e", "_west" => "w",
    "chimp" => "C", "neanderthal" => "N", "denisovan" => "D")
end
# switch some minor/major partners
setgamma!(altg[4].edge[21], 0.501)
setgamma!(altg[2].edge[23], 0.501)
setgamma!(altg[1].edge[23], 0.501)
setgamma!(altg[1].edge[3],  0.501)
for i in [-4,-5,-10,-14,-16] rotate!(altg[4], i); end
for i in [-4,-6,-9,-10,-11,-14] rotate!(altg[3], i); end
for i in [-4,-9,-8,-11,-12,-13] rotate!(altg[2], i); end
for i in [-3,-4,-5,-6,-13] rotate!(altg[1], i); end
for i in [-2,-4,-5,-10,-12] rotate!(altg[5], i); end

pargs = (tipcex=0.8, tipoffset=0.1, style=:majortree, arrowlen=0.1, majorhybridedgecolor = "deepskyblue")
R"cairo_pdf"("figures/SMfig_nets_samefitasg4.pdf", height=5, width=9);
R"layout"([1 2 3 4; 1 5 6 7]);
R"par"(mar=[0,0,0,0]);
plot(fleg_nobott; pargs...);
R"mtext"(text="g4", line=-3, adj=0.1)
plot(g4_l4_c34_nobot; pargs...);
R"mtext"(text="g4-l4-c3444", line=-3, adj=0.1, cex=0.6)
for (net,irep) in zip(altg, altr)
  plot(net; pargs...)
  R"mtext"(text="replicate $irep", line=-3, adj=0.1, cex=0.6)
end
R"dev.off"();

end # of `if false`
