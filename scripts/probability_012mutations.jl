#= Calculation of the marginal probability of 0, 1, 2+ mutations per site,
   across the entire gene tree.

This gene tree is random from the coalescent, so these probabilities are
calculated from each simulated gene tree, then averaged across gene trees.

Given each gene tree, this probability depends on the site's rate, assuming a
distribution of rates across sites that is a discretized-Gamma with 10 rate
categories. The shape of the Gamma distribution is α=0.35 by default, but
functions can use α as an input parameter.

example use, easy (sequentially) but somewhat time-consuming:

```julia
include("../scripts/probability_012mutations.jl")
df = probability012more()
CSV.write("probability012more_byrep.csv", df)
```

then summarize, by combining replicates from each parameter combination:

```julia
df = probability012more_summarize("probability012more_byrep.csv")
CSV.write("probability012more_byparam.csv", df)
```

plot these data:

```julia
using CSV, DataFrames, AlgebraOfGraphics, CairoMakie
p0 = data(df) * mapping(
    :nind => nonnumeric => "number of individuals",
    :p2more_cond_mean => "P{2+ mutations | 1+ mutation}",
    marker = :mutrate => nonnumeric,
    color = :lindist => nonnumeric => "rate variation across lineages",
    dodge_x = :net,
    ) * visual(Scatter, alpha=0.8)
draw(p0, scales(DodgeX=(; width=0.5), Color=(; palette=from_continuous(:viridis))))
```

conclusion: factors causing P{2 or more mutations | 1+ mutation} to vary are:
  + mutation rate
  + number of individuals

```julia
dfsum = sort!(combine(groupby(df, [:nind, :mutrate]),
  :p0_mean     => mean => :p0,
  :p1_mean     => mean => :p1,
  :p2more_mean => mean => :p2more,
  :p2more_cond_mean => mean => :p2more_cond
  ), [:mutrate, :nind])
```

 Row │ nind   mutrate  p0        p1         p2more       p2more_cond 
     │ Int64  Float64  Float64   Float64    Float64      Float64     
─────┼───────────────────────────────────────────────────────────────
   1 │     1  1.25e-8  0.989466  0.0103397  0.000194575    0.0177393
   2 │     2  1.25e-8  0.987432  0.0122941  0.000274029    0.0210915
   3 │    10  1.25e-8  0.981096  0.0182942  0.00061019     0.0315321
   4 │     1  1.25e-7  0.909293  0.076267   0.0144396      0.15538
   5 │     2  1.25e-7  0.89591   0.0852163  0.0188741      0.178056
   6 │    10  1.25e-7  0.855307  0.108562   0.0361303      0.246897

=#

using PhyloNetworks
using Distributions
using DataFrames, CSV

"""
    hasgenetrees_substitutions(dir::AbstractString)

`true` if `dir` contains a file named "001/gts_nsubst": the file with gene trees
used as input for `seqgen`, with branch lengths in substitutions per site.
`false` otherwise.
"""
function hasgenetrees_substitutions(dir::AbstractString)
    file = joinpath(dir,"001","gts_nsubst")
    return (isfile(file) && filesize(file)>0)
end

"""
    parse_folder2parameters(dir::AbstractString)

Named tuple of parameters that `dir` contains output for.

# example

```julia
julia> tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate";

julia> parse_folder2parameters(tmp)
(net = "fleg-pruned", nind = 1, nbiallelic = 10000, lindist = 0.0, mutrate = 1.25e-7)
```
"""
function parse_folder2parameters(dir::AbstractString)
    paramstr = split(dir, "_")
    net = paramstr[1] # fleg or fleg-pruned
    paramdict = Dict(x[2] => x[1] for x in rsplit.(paramstr[2:end],"-",limit=2))
    nind    = parse(Int,     paramdict["indls"])
    nbiallelic = parse(Int,  paramdict["biall"])
    lindist = parse(Float64, paramdict["lindist"])
    mutrate = parse(Float64, paramdict["subrate"])
    return (net=net, nind=nind, nbiallelic=nbiallelic, lindist=lindist, mutrate=mutrate)
end

"""
    readgenetrees_seqgenformat(file::AbstractString)

Vector of gene trees in `file`, formatted for seqgen with 1 tree per line
and the number of desired sites to be simulated before the tree like this:
[500]((a:1,b:1):1,c:2)

# example

```julia
tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate/001/gts_nsubst"
readgenetrees_seqgenformat(tmp)
```
"""
function readgenetrees_seqgenformat(file::AbstractString)
    vnet = HybridNetwork[]
    open(file) do s
        for line in eachline(s)
            line = replace(line, r"^\[\d+\]" => s"") # removes [500]
            isempty(line) && continue
            push!(vnet, readTopology(line,false)) # verbose=false
        end
    end
    return vnet
end

totaledgelength(net) = sum(e.length for e in net.edge)

"""
    discretegammarates(α=0.35)

Vector of 10 rate categories for the discrete Gamma distribution of rates
across sites when sequences are simulated with `seqgen` with option
`-g 10` to ask for Gamma-distributed rates and 10 rate categories.
The shape parameter α was used for seqgen, with option `-a 0.35`.

The gamma distribution is discretized by taking 10 intervals of probability 1/10
each, then taking the mean within each interval. This is as done by seqgen.
To obtain the intervals of equal probabilities, obtain the integral under each
interval (to get the interval means), we use the `Gamma` distribution from
[Distributions](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Gamma).

The mean of the censored Gamma is not implemented, but we can use the following
trick. Let X₀ ~ Γ(α,1) and X₁ ~ Γ(α+1,1). Then, using Γ(α+1) = αΓ(α) and the
density of X₀ (exp(-t) t^{α-1} dt), we get:
E(X₀ indicator{X₀ < x}) = α P{X₁ < αx}

For the rescaled Y₀ = X₀/α. We get: E(Y₀ indicator{Y₀ < y}) = P{X₁ < αy}.
This leads to the conditional mean of Y₀~Γ(α,α) based on the cdf of X₁~Γ(α+1,1):

    E(Y₀ | y<Y₀<z) = (P{X₁ < αz} - P{X₁ < αy}) / P{y<Y₀<z}

# example

```julia
julia> discretegammarates()
10-element Vector{Float64}:
 0.0007408197161752766
 0.010023497386830775
 0.04101565373067346
 0.10733088085014231
 0.22530347862560146
 0.4182656216521762
 0.7256397696095872
 1.2278652063396798
 2.1443594991397714
 5.099455572949363
```
"""
function discretegammarates(α=0.35)
    dY = Gamma(α, 1/α) # scale θ = 1/α
    dX = Gamma(α+1, 1)
    qY = quantile.(dY, 0:0.1:1) # 11 values: includes 0 and Inf
    # [0.0028583487040634093,0.020807468876047546,0.0670688315915516,0.15609529648641018,0.3067699596052665,0.5482525405438332,0.9338020936265511,1.5833999842821815,2.8877008177391517]
    pX = cdf.(dX, α .* qY)
    rates = (pX[2:11] - pX[1:10]) .* 10
    all(qY[1:10] .< rates .< qY[2:11]) ||
        error("rate categories are not within their interval bounds")
    return rates
end

const discretegammarates_α035 = [
 0.0007408197161752766, 0.010023497386830775, 0.04101565373067346,
 0.10733088085014231, 0.22530347862560146, 0.4182656216521762,
 0.7256397696095872, 1.2278652063396798, 2.1443594991397714, 5.099455572949363]

"""
    probability012more(treelength::Number, rates::AbstractVector)

Named tuple of P(0 mutations), P(1 mutation) and P(2+ mutations)
for a gene tree of total length `treelength` and under a model of rate
variation with 10 `rates`, each with probability 1/10.

`rates` should be a vector of size 10. Each of the 10 rates is assumed to have
probability 1/10.

# example

```julia
julia> probability012more(0.5, discretegammarates_α035)
(p0 = 0.7284124447705469, p1 = 0.14974235464990363, p2more = 0.12184520057954944)
```
"""
function probability012more(treelength::Real, rates::AbstractVector)
    glen = treelength .* rates
    # probabilities conditional on the mutation rate
    p0 = exp.(-glen)    # 0 mutations: Poisson(0) = exp(-λ)
    p1 = p0 .* glen    # 1 mutation:  Poisson(1) = λexp(-λ)
    p2 = 1 .- p0 .- p1 # 2 or more mutations
    # marginal probabilities
    return (p0 = sum(p0)/10, p1 = sum(p1)/10, p2more = sum(p2)/10)
end

"""
    probability012more(treefile::AbstractString, rates::AbstractVector)

Named tuple of P(0 mutations), P(1 mutation) and P(2+ mutations), where
probabilities are averaged over gene trees found in file ``
of directory `dir. In this gene tree file:
- branch lengths are assumed in **substitutions** per site
- the format is for `seqgen` with a number of desired site in front of each tree,
  see `readgenetrees_seqgenformat`.
- each tree is given equal weight, that is, the same number of sites is assumed
  to come from each tree.
`
`rates` should be a vector of size 10. Each of the 10 rates is assumed to have
probability 1/10.

# example

```julia
julia> tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate/001/gts_nsubst";

julia> probability012more(tmp, discretegammarates_α035)
(p0 = 0.9095664768614613, p1 = 0.07647390198461576, p2more = 0.013959621153926441)
```
"""
function probability012more(treefile::AbstractString, rates::AbstractVector)
    genetree = readgenetrees_seqgenformat(treefile)
    p0,p1,p2 = 0.0,0.0,0.0
    for gt in genetree
        p = probability012more(totaledgelength(gt), rates)
        p0 += p[:p0]
        p1 += p[:p1]
        p2 += p[:p2more]
    end
    nt = length(genetree)
    return (p0 = p0/nt, p1 = p1/nt, p2more = p2/nt)
end

"""
    probability012more!(df::DataFrame, dir::AbstractString, rates::AbstractVector)

Add to data frame `df` one row for each replicate data set in directory `dir`,
containing the parameter set for the directory, replicate number, and
marginal probabilities of 0, 1, and 2 or more mutations per site.

# example

```julia
julia> tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate";

julia> tmp_df = DataFrame(net=String[], # initialize an empty data frame
        nind=Int[], nbiallelic=Int[], lindist=Float64[], mutrate=Float64[],
        irep=Int[], p0 = Float64[], p1 = Float64[], p2more = Float64[]);

julia> probability012more!(tmp_df, tmp, discretegammarates_α035)
[ Info: 1.370 mins for fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate
100×9 DataFrame
 Row │ net          nind   nbiallelic  lindist  mutrate  irep   p0        p1         p2more    
     │ String       Int64  Int64       Float64  Float64  Int64  Float64   Float64    Float64   
─────┼─────────────────────────────────────────────────────────────────────────────────────────
   1 │ fleg-pruned      1       10000      0.0  1.25e-7      1  0.909566  0.0764739  0.0139596
   2 │ fleg-pruned      1       10000      0.0  1.25e-7      2  0.909595  0.0764499  0.0139548
  ⋮  │      ⋮         ⋮        ⋮          ⋮        ⋮       ⋮       ⋮          ⋮          ⋮
  99 │ fleg-pruned      1       10000      0.0  1.25e-7     99  0.909372  0.0766054  0.014023
 100 │ fleg-pruned      1       10000      0.0  1.25e-7    100  0.909517  0.0765032  0.0139797
```
"""
function probability012more!(df::DataFrame, dir::AbstractString, rates::AbstractVector)
    net, nind, nbiallelic, lindist, mutrate = parse_folder2parameters(dir)
    repdir = filter(contains(r"\d+"), readdir(dir))
    starttime = time()
    for repstr in repdir
        irep = parse(Int, repstr)
        treefile = joinpath(dir, repstr, "gts_nsubst")
        probs = probability012more(treefile, rates)
        push!(df, (net,nind,nbiallelic,lindist,mutrate, irep,probs...))
    end
    endtime = time()
    elapsed = round((endtime - starttime) / 60, digits=3)
    @info "$elapsed mins for $dir"
    GC.gc()
    return df
end

"""
    probability012more(
        paramfilters::Vector{Regex}=Regex[],
        rates::AbstractVector=discretegammarates_α035
    )

Data frame with one row for each replicate data set from each parameter
combination in the current directory, whose folder name matches each
regular expression in `paramfilters`.
If a parameter combination folder does not contains a non-empty file
"001/gts_nsubst", then this folder is skipped.

Columns: parameters, replicate number, and marginal probabilities of 0, 1,
and 2 or more mutations per site.

# example

```julia
julia> paramf = ["fleg-pruned_", "_1.25e-7-subrate"];

julia> df = probability012more(Regex.(paramf))
[ Info: starting: 30 parameter combination folders
[ Info: 1.48 mins for fleg_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-8-subrate
...
3000×9 DataFrame
  Row │ net          nind   nbiallelic  lindist  mutrate  irep   p0        p1         p2more     
      │ String       Int64  Int64       Float64  Float64  Int64  Float64   Float64    Float64    
──────┼──────────────────────────────────────────────────────────────────────────────────────────
    1 │ fleg-pruned      1       10000      0.0  1.25e-7      1  0.909566  0.0764739  0.0139596
  ⋮   │      ⋮         ⋮        ⋮          ⋮        ⋮       ⋮       ⋮          ⋮          ⋮

julia> filename = "probability012more_" * join(replace.(paramf, "_" => ""), '_') * ".csv"
"probability012more_fleg-pruned_1.25e-7-subrate.csv"

julia> CSV.write(joinpath("/nobackup2/ane/", filename), df);
```
"""
function probability012more(
    paramfilters::Vector{Regex}=Regex[],
    rates::AbstractVector=discretegammarates_α035
)
    df = DataFrame(net=String[],
            nind=Int[], nbiallelic=Int[], lindist=Float64[], mutrate=Float64[],
            irep=Int[], p0 = Float64[], p1 = Float64[], p2more = Float64[])
    paramset_dirs = filter(contains(r"\w+_\d.+"), readdir())
    for f in paramfilters
        filter!(contains(f), paramset_dirs)
    end
    @info "starting: $(length(paramset_dirs)) parameter combination folders"
    for dir in paramset_dirs
        if !hasgenetrees_substitutions(dir)
            @warn "no 'gts_nsubst' file: skip $dir"
            continue
        end
        probability012more!(df, dir, rates)
    end
    return df
end


using Statistics

"""
    probability012more_summarize(df::DataFrame)
    probability012more_summarize(csv_filename::AbstractString)

Average probabilities of 0,1,2+ mutations/site over replicates, also with
the conditional probability:

P{2+ mutations | 1+ mutation} = P{2+} / (P{1 mutation} + P{2+ mutations})

This is an overestimate of the conditional probability:

P{2+ mutations | variable site}

because it is possible that a site with 2+ mutations looks invariable at the
leaves: 1 mutation can be reversed in all its descendants by one (or multiple
separate) back-mutation(s). However, a site with 0 mutations is necessarily
invariable: P{variable site} > P{1+ mutations}

# example

```julia
julia> res = probability012more_summarize("probability012more_byrep.csv");

julia> CSV.write("probability012more_byparam.csv", res);

julia> sort!(res, [:p2more_cond_mean]; rev=[true]);
```
"""
function probability012more_summarize(filename::AbstractString)
    df = CSV.read(filename, DataFrame)
    probability012more_summarize(df)
end
function probability012more_summarize(df::DataFrame)
    # add column for the conditional probability: P{2+ mutations | 1+ mutation}
    df.p2more_cond = df.p2more./(df.p2more .+ df.p1)
    # do *not* group by 'nbiallelic' bc probabilities do not depend on # of sites
    gdf = groupby(df, [:net, :nind, :lindist, :mutrate])
    @info "$(length(keys(gdf))) parameter combinations"
    tmp = combine(gdf, nrow => "nreps")
    if !all(tmp[!,:nreps] .== 200) # 100 * 2 conditions for nbiallelic
        @warn "some parameter combinations do not have 100 reps"
        subset!(tmp, :nreps => x -> x .!== 200)
        @show tmp
    end
    res = combine(gdf, :p0 => mean, :p0 => std, :p1 => mean, :p1 => std,
      :p2more => mean, :p2more => std, :p2more_cond => mean, :p2more_cond => std)
    return res
end
