# Software installation notes

Organized by software package.

## PhyloNetworks - done on Franklin

Julia 1.10.0 available on darwins/franklins as of 2024-2-6.

note: might be able to disregard #1 below, as `Manifest.toml` version documentation has improved in recent major updates. 

1. Set up Franklin Julia environment by transferring local Manifest.toml (from correct version) to `~/.julia/environments/v1.10/`
2. Activate julia from `f4-ratevar` directory with `julia --project` to download packages from Project.toml like so:

```{julia package mode}
activate .
instantiate
```

To get latest commit (not included in latest tagged version v0.16.3):
(commit is dated 2024-1-5)

```{julia package mode}
add PhyloNetworks#f46226c
```

## R packages

R version 4.4.0 (2024-4-24) -- "Puppy Cup"
installed on the cluster, default in `/usr/bin/R` as of 2024-5-20.

Install packages within R, locally to `~/R/x86_64-pc-linux-gnu-library/4.4`:
** Note that when major R updates occur, packages won't be found, and need to be re-installed (or moved?).

### SiPhyNetwork - done on Franklin

version 1.1.0 downloaded as of 2024-5-20.

```{r}
install.packages("devtools") # takes forever!!
install.packages('ape')
devtools::install_github("jjustison/SiPhyNetwork")
```

### admixtools - done on Franklin

version 2.0.4 git commit a3d1d5b, downloaded as of 2024-5-20.

```{r}
devtools::install_github("uqrmaie1/admixtools")
```

### updating

```{r}
library(remotes)
remotes::update_packages()
```
