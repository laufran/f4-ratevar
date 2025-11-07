# Software installation notes

Organized by software package.

## PhyloNetworks - done on Franklin

Julia 1.11.1 available on darwins/franklins as of 2025-5-14.

Activate julia from `f4-ratevar` directory with `julia --project` to download packages from Project.toml like so:

```{julia package mode}
activate .
instantiate
```

## R packages

R version 4.5.0 (2025-4-24) -- "Puppy Cup"
installed on the cluster, default in `/usr/bin/R` as of 2025-4-19.

Install packages within R, locally to `~/R/x86_64-pc-linux-gnu-library/4.5`:
** Note that when major R updates occur, packages won't be found, and need to be re-installed (or moved?).

### admixtools - done on Franklin

version 2.0.4 git commit a3d1d5b, downloaded as of 2024-5-20. (for g1/g4)
version 2.0.9 git commit d394bd0 (for additional graphs beyond g1/g4)

```{r}
install.packages("devtools") # takes forever!!
devtools::install_github("uqrmaie1/admixtools")
```

### updating

```{r}
library(remotes)
remotes::update_packages()
```
