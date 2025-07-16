# pipeline to study lineage rate variation's effect on f4 stats

Code to reproduce the simulation study in: Frankel & Ané (2025) "Low accuracy of complex admixture graph inference from f-statistics".

Preprint available here: <https://www.biorxiv.org/content/10.1101/2025.03.07.642126v2>

## goal and model

Goal: simulate sequence data (SNPs) under a network model reflecting
archaic human population genetics with ILS and lineage rate variation,
to answer these questions:

1. How does f4 perform when there is lineage rate variation?
2. How much lineage rate variation causes consequential type-1 error,
or a reduction in power?
3. How are the admixture graphs inferred from f4 results affected?

## dependencies

- julia packages: see `Project.toml`
- seqgen for the simulation of sequences: Seq-Gen v1.3.4 downloaded from
  [github](git@github.com:rambaut/Seq-Gen.git) at [commit from 2019-08-29](https://github.com/rambaut/Seq-Gen/commit/e8660d73769297d25a88b40d5b828a6fce7784b5).
- R packages:
  - [admixtools](https://github.com/uqrmaie1/admixtools), see [notes/software_install.md](notes/software_install.md#admixtools---done-on-franklin)
