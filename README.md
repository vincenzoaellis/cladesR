# cladesR

<!-- badges: start -->
<!-- badges: end -->

> R functions for identifying clades in a phylogenetic tree that are larger or smaller than expected based on a homogeneous birth-death model.

## Overview

cladesR provides R functions for finding clades in a dated phylogeny that are
larger or smaller than expected under a homogeneous birth–death model of
diversification. This package allows for reproduction of analyses in Nee et al.
(1992) and Ricklefs (2003, 2014, etc.). You will need an ultrametric phylogeny,
then you can slice it at a given time point in the past, and finally examine the
resulting clades. There are also some plotting functions to keep track of the
clade assignments of lineages. I've written this with Claude Code (Opus 4.8), so
perhaps extra caution is required before using; although I've done some pencil
and paper examples and the functions seem to work as intended. If you notice any
problems or can think of improvements, please open an issue.

## Installation

You can install the development version of cladesR from
[GitHub](https://github.com/vincenzoaellis/cladesR) with:

``` r
# install.packages("remotes")
remotes::install_github("vincenzoaellis/cladesR")
```

cladesR depends on **ape** and **phytools**; **ggplot2** is optional and only
needed for the plotting functions.

## Functions

| Function | What it does |
|----------|--------------|
| `clade_tests()` | Per-clade one-sided tests vs. the geometric null (`too_few` / `too_many`), with optional multiple-testing correction |
| `sd_sim_test()` | Tests the observed SD of clade sizes against the broken-stick (conditional) null of Nee et al. 1992 |
| `geom_expectation()` | Summary statistics + theoretical geometric expectations (SD, SD/mean, Gini–Simpson, proportion singletons) |
| `clade_rank_data()` | ln(rank)-vs-size data and the geometric prediction line (optional ggplot) |
| `extract_clade_sizes()` | Slice an ultrametric tree at a given age and return the resulting clades, sizes, species, and nodes |
| `build_clade_table()` | Species × time-slice clade-membership table (ready for nested ANOVA) |
| `ancestor_descendant_stats()` | Parent→descendant clade-size statistics across adjacent slices |
| `compute_pulse_score()` | Per-tip index of how consistently a lineage sits in over- (or under-) sized clades across many slices |
| `plot_cladetracker()` | Multi-panel "clade tracker" plot linking collapsed clades across time slices |

## Example

### A clade that fits the model (simulated)

A pure-birth (Yule) tree is the textbook case: equal-aged clades follow a
geometric distribution by construction (Kendall 1948), so the diagnostics should
*fail to reject* the null.

``` r
library(cladesR)

# An ultrametric, 30-"Ma"-tall pure-birth tree.
set.seed(17)
tree <- phytools::pbtree(b = 1, n = 150, scale = 30)

# Slice 10 Ma before the present and collect the surviving clades.
slice <- extract_clade_sizes(tree, age_before_present = 10)

# Summary statistics and the theoretical geometric expectations.
geom_expectation(slice$clade_sizes)

# Per-clade tests vs the geometric null -- nothing is flagged here.
clade_tests(slice$clade_sizes)

# Is the spread of clade sizes consistent with the null? (broken-stick test)
sd_sim_test(slice$clade_sizes)

# ln(rank)-vs-size, with the geometric prediction line.
clade_rank_data(slice$clade_sizes, plot = TRUE)
```

### A clade that does not (Furnariidae, via `clootl`)

Real radiations often depart from the homogeneous model. Here we pull the
ovenbirds and woodcreepers (the group analysed by Ricklefs 2014) from the
[`clootl`](https://github.com/eliotmiller/clootl) bird phylogeny, which is dated
in millions of years.

``` r
library(clootl)

taxa <- taxonomyGet(taxonomy_year = 2025)
furn <- taxa[taxa$FAMILY == "Furnariidae (Ovenbirds and Woodcreepers)", ]
phy  <- extractTree(species = furn$sci_name_2025)

# At a single 5-Ma slice, several clades are larger than the geometric null
# predicts -- the signature of diversification pulses.
furn5 <- extract_clade_sizes(phy, age_before_present = 5)
subset(clade_tests(furn5$clade_sizes), too_many)

# Multi-slice workflow: slice at several ages, track how clade sizes change
# from ancestors to descendants, and score every tip for how often it sits in
# an oversized clade.
ages   <- c(15, 10, 5)
slices <- lapply(ages, function(a) extract_clade_sizes(phy, age_before_present = a))

ancestor_descendant_stats(build_clade_table(slices, ages))
compute_pulse_score(phy, slice_ages = ages)

# Visualise the clade assignments through time.
plot_cladetracker(phy, slice_list = slices, slice_ages = ages)
```

### Partitioning trait variance across clade levels

`build_clade_table()` turns the slices into a taxonomy-style table: one row per
species, one column per time slice, each cell holding the monophyletic clade the
species belongs to at that depth. Because the columns are nested grouping
factors, you can partition the variance of a trait across phylogenetic levels --
a clade-based, equal-aged alternative to a nested ANOVA on Linnaean ranks
(Ricklefs 2005). Clades are a random sample of possible clades, so the
statistically appropriate model is a **mixed model with nested random effects**,
which yields the variance component at each level directly.

``` r
ctab <- build_clade_table(slices, ages)
head(ctab)
#>                   species clade_15Ma clade_10Ma clade_5Ma
#> 1 Cranioleuca albicapilla   clade_01   clade_01 clade_001
#> 2    Cranioleuca albiceps   clade_01   clade_01 clade_001
#> 3 Cranioleuca antisiensis   clade_01   clade_01 clade_001
#> ...

# A simulated, normally distributed trait (substitute a real measurement here).
set.seed(42)
ctab$trait <- rnorm(nrow(ctab))

# Nested random effects -> variance component (and % of total) at each depth.
library(lme4)
fit <- lmer(trait ~ 1 + (1 | clade_15Ma / clade_10Ma / clade_5Ma), data = ctab)
vc  <- as.data.frame(VarCorr(fit))
data.frame(level = vc$grp, variance = round(vc$vcov, 4),
           pct = round(100 * vc$vcov / sum(vc$vcov), 1))
#>                               level variance pct
#> 1 clade_5Ma:(clade_10Ma:clade_15Ma)   0.0000   0
#> 2             clade_10Ma:clade_15Ma   0.0000   0
#> 3                        clade_15Ma   0.0000   0
#> 4                          Residual   0.9416 100
```

Because this trait is just random noise with no phylogenetic structure, **all of
the variance is residual** -- here the residual is the variance *among lineages
that share the same base-level (youngest, 5 Ma) clade*, i.e. variation not
explained by clade membership at any slice. (lme4 reports a "singular fit", as
expected when the clade-level variances are estimated at zero.) A real,
phylogenetically structured trait would instead push variance up into the
deeper, older clade levels. With no extra dependencies the same decomposition is
available in base R via error strata:
`aov(trait ~ Error(clade_15Ma / clade_10Ma / clade_5Ma), data = ctab)`.

## References

- Nee S, Mooers AO, Harvey PH (1992) Tempo and mode of evolution revealed from
  molecular phylogenies. *PNAS* 89:8322–8326.
  <https://doi.org/10.1073/pnas.89.17.8322>
- Ricklefs RE (2003) Global diversification rates of passerine birds. *Proc. R.
  Soc. B* 270:2285–2291. <https://doi.org/10.1098/rspb.2003.2489>
- Ricklefs RE (2005) Small clades at the periphery of passerine morphological
  space. *Am. Nat.* 165:651–659. <https://doi.org/10.1086/429676>
- Ricklefs RE (2014) Reconciling diversification: random pulse models of
  speciation and extinction. *Am. Nat.* 184:268–276.
  <https://doi.org/10.1086/676642>
- Kendall DG (1948) On the generalized "birth-and-death" process. *Ann. Math.
  Stat.* 19:1–15. <https://doi.org/10.1214/aoms/1177730285>

## License

GPL-2.
