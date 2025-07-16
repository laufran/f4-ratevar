#= Estimation of the probabilities
   P{site is variable} and P{site is multiallelic | variable}

using the number of sites that were simulated to get 1 variable site per gene,
and number of genes to be simulated (1 variable site/gene) to get the desired
(10,000 or 100,000) biallelic sites total.

in the folder for each replicate:
- `seqgen_log.txt`: csv file with a header and 1 line.
  `sim_multi` = number of multiallelic sites saved to vcf, until the desired
   number of biallelic sites were simulated (10_000 or 100_000)
- `gene_log.txt`:  csv file with 2 columns `gene_number,ninvar_sites`
   and 1 line per gene, for the genes that were "used":
   `100_000` + `sim_multi` genes.
  `ninvar_sites` = number of sites until 1 variable site was generated.

example use, quite fast:

```julia
include("../scripts/probability_multiallelic.jl")
df = probability_variable_multiallelic()
CSV.write("probability_multiallelic.csv", df)
```

plot these data:

```julia
using CSV, DataFrames, AlgebraOfGraphics, CairoMakie
pv = data(df) * mapping(
    :nind => nonnumeric => "number of individuals",
    :pvar => "P{variable}",
    marker = :mutrate => nonnumeric,
    color = :lindist => nonnumeric => "rate variation across lineages",
    dodge_x = :net) * visual(Scatter, alpha=0.8)
pm = data(df) * mapping(
    :nind => nonnumeric => "number of individuals",
    :pmulti_mean => "P{multiallelic | variable}",
    marker = :mutrate => nonnumeric,
    color = :lindist => nonnumeric => "rate variation across lineages",
    dodge_x = :net) * visual(Scatter, alpha=0.8)
fig = draw(pv, scales(DodgeX=(; width=0.5),
  Color=(; palette=from_continuous(:viridis))));
fig.figure.layout # GridLayout[1:1, 1:2]
draw!(fig.figure[2,1], pm, scales(DodgeX=(; width=0.5),
  Color=(; palette=from_continuous(:viridis))));
legend = content(fig.figure[1,2])
fig.figure[:,2] = legend
fig
save("probabilitymultiallelic.pdf", fig)
```

conclusion: factors causing P{multiallelic site | variable site} to vary are:
  + mutation rate
  + number of individuals

```julia
dfsum = sort!(combine(groupby(df, [:nind, :mutrate]),
  :pvar => mean => :pvar,
  :pmulti      => mean => :pmulti,
  :pmulti_mean => mean => :pmulti_mean,
  ), [:mutrate, :nind])
```

 Row │ nind   mutrate  pvar       pmulti      pmulti_mean 
     │ Int64  Float64  Float64    Float64     Float64     
─────┼────────────────────────────────────────────────────
   1 │     1  1.25e-8  0.0104916  0.00503092   0.00503006
   2 │     2  1.25e-8  0.0125069  0.0067808    0.00677933
   3 │    10  1.25e-8  0.0188486  0.0117114    0.0117081
   4 │     1  1.25e-7  0.0881653  0.0473741    0.0473318
   5 │     2  1.25e-7  0.101898   0.0621713    0.0621055
   6 │    10  1.25e-7  0.142798   0.102446     0.102308

=#

# load script that defines parse_folder2parameters:
include("probability_012mutations.jl")

"""
    probability_variable(folder)

Probability that a site is variable, estimated from a csv file in `folder` with
2 columns `gene_number,ninvar_sites` and 1 line per gene, where `ninvar_sites`
is the number of sites that had to be simulated to get 1 variable site.
We want the marginal probability P{variable}, averaged across all genes.

Conditionally on the gene tree (with branch lengths in substitutions/site),
`ninvar_sites` has a geometric distribution. But it does *not* marginally,
when we don't condition on the gene tree. So we use the partial information:
P(ninvar_sites = 1 | gene) = P(variable | gene)
because this relationship is maintained if we integrate out the gene tree:

P(variable) is estimated as the proportion of gene trees for which we observe
`ninvar_sites` = 1

# example

```julia
julia> tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate/001";

julia> probability_variable(tmp)
(p_variable = 0.08914839077185988, n_firstsite = 939, n_trials = 10533)
```
"""
function probability_variable(dir)
    file = joinpath(dir, "gene_log.txt")
    ns = CSV.read(file, NamedTuple)[:ninvar_sites]
    ngenes = length(ns)
    filter!(isequal(1), ns)
    return (p_variable=length(ns)/ngenes, n_firstsite=length(ns), n_trials=ngenes)
end
"""
    probability_multiallelic(folder)

Probability that a site is multiallelic (3+ observed states) conditional on this
site being variable, estimated from the number of genes that had to be simulated
to get to a desired number of biallelic sites. From each gene, the first variable
site was saved. This site is multiallelic with the probability of interest:
p = P(multiallelic | variable).
The number of genes N to simulate to get 10_000 biallelic sites has a shifted
negative binomial distribution: N - 10_000 ~ negbinom(r=10_000, 1-p)
whose mean is μ = 10_000 * p/(1-p), variance μ/(1-p).
The maximum likelihood estimate is from the method of moments:
μ̂ ≈ number of genes that gave a multiallelic site, then invert:
p̂ = μ̂/(μ̂ + 10_000) --- but this is biased. A minimum-variance unbiased estimate is:

p̂ = μ̂/(μ̂ + 10_000-1)

# example

```julia
julia> tmp = "fleg-pruned_1-indls_11800-genes_10000-biall_0.0-lindist_1.25e-7-subrate/001";

julia> probability_multiallelic(tmp)
(p_multiallelic = 0.05060767185719711, n_total = 10533, n_multi = 533, enough = true)
```
"""
function probability_multiallelic(dir)
    file = joinpath(dir, "seqgen_log.txt")
    # @info "  starting: $dir"
    #= some files are corrupt, starting with lots of junk. Only the last 2 lines are correct.
    nt = CSV.read(file, NamedTuple)
    # (input_genes = [11800], req_bi = [10000], sim_bi = [10000], sim_multi = [533], failed_genes = [0])
    nreq = nt[:req_bi][1]
    nbi  = nt[:sim_bi][1]
    nmul = nt[:sim_multi][1]
    =#
    header = "input_genes,req_bi,sim_bi,sim_multi,failed_genes"
    dataline = ""
    readnextline = false
    for line in eachline(file)
        strip(line)
        if readnextline
            dataline = line
        end
        readnextline = line == header
    end
    data = split(dataline, ",")
    nreq = parse(Int, data[2])
    nbi  = parse(Int, data[3])
    nmul = parse(Int, data[4])
    ntot = nbi + nmul
    denom = ntot-1
    enough = nreq == nbi
    if !enough
        @warn "there weren't enough gene trees: will remove the bias correction"
        denom = ntot
    end
    pm = nmul / denom
    return (p_multiallelic=pm, n_total=ntot, n_multi=nmul, enough=enough)
end

function probability_variable_multiallelic!(df::DataFrame, dir::AbstractString)
    net, nind, nbiallelic, lindist, mutrate = parse_folder2parameters(dir)
    repdir = filter(contains(r"\d+"), readdir(dir))
    @info "starting: $dir"
    # starttime = time()
    p_variable = Float64[]; n_firstsite = Int[]; n_trials = Int[]
    p_multi = Float64[]; ngenes_sim = Int[]; ngenes_multi = Int[]; enough = Bool[]
    for repstr in repdir
        # irep = parse(Int, repstr)
        repdir = joinpath(dir, repstr)
        res = probability_variable(repdir)
        push!(p_variable,  res[:p_variable])
        push!(n_firstsite, res[:n_firstsite])
        push!(n_trials,    res[:n_trials])
        res = probability_multiallelic(repdir)
        push!(p_multi, res[:p_multiallelic])
        push!(ngenes_sim,  res[:n_total])
        push!(ngenes_multi,res[:n_multi])
        push!(enough,      res[:enough])
    end
    pv = sum(n_firstsite)/sum(n_trials)
    n_multi_tot = sum(ngenes_multi)
    n_enough = sum(enough)
    n_not_enough = sum(.!enough)
    n_not_enough>0 && @info "not enough genes in $n_not_enough reps, in $dir"
    pm1 = n_multi_tot / (sum(ngenes_sim) - n_enough)
    pm2 = mean(p_multi)
    push!(df, (net,nind,nbiallelic,lindist,mutrate, pv, pm1,pm2))
    # endtime = time()
    # elapsed = round((endtime - starttime) / 60, digits=3)
    # @info "$elapsed mins for $dir"
    GC.gc()
    return df
end

"""

# example

```julia
julia> paramf = ["fleg-pruned_", "_1.25e-7-subrate", "100000-biall"];

julia> df = probability_variable_multiallelic(Regex.(paramf))
[ Info: starting: 15 parameter combination folders
[ Info: starting: fleg-pruned_1-indls_118000-genes_100000-biall_0.0-lindist_1.25e-7-subrate
...

julia> filename = "probability_multiallelic_" * join(replace.(paramf, "_" => ""), '_') * ".csv";

julia> CSV.write(joinpath("/nobackup2/ane/", filename), df);
```
"""
function probability_variable_multiallelic(
    paramfilters::Vector{Regex}=Regex[]
)
    df = DataFrame(net=String[],
            nind=Int[], nbiallelic=Int[], lindist=Float64[], mutrate=Float64[],
            pvar=Float64[], pmulti=Float64[], pmulti_mean=Float64[])
    paramset_dirs = filter(contains(r"\w+_\d.+"), readdir())
    for f in paramfilters
        filter!(contains(f), paramset_dirs)
    end
    @info "starting: $(length(paramset_dirs)) parameter combination folders"
    for dir in paramset_dirs
        probability_variable_multiallelic!(df, dir)
    end
    return df
end
