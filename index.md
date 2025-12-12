# ggrecipes

A collection of easy-to-use custom
[ggplot2](https://ggplot2.tidyverse.org/index.html)-based functions for
data exploration and analysis. Each function handles data preprocessing
and returns a ggplot2 object that can be further customized using
standard ggplot2 syntax. Includes general-purpose and domain-specific
visualizations.

## Features

![](reference/figures/demo.gif)

**General:**

General-purpose visualizations

- [`gg_splitcorr()`](https://ignophi.github.io/ggrecipes/reference/gg_splitcorr.html) -
  Split-Correlation Heatmap
- [`gg_rankshift()`](https://ignophi.github.io/ggrecipes/reference/gg_rankshift.html) -
  Paired box-/bar- plots with rank change
- [`gg_criteria()`](https://ignophi.github.io/ggrecipes/reference/gg_criteria.html) -
  Criteria heatmap with optional barplots
- [`gg_conf()`](https://ignophi.github.io/ggrecipes/reference/gg_conf.html) -
  Confusion/contingency table bubble plot

**Bioinformatics:**

Sequence analysis and genomics visualizations

- [`gg_geno()`](https://ignophi.github.io/ggrecipes/reference/gg_geno.html) -
  Biallelic genotype visualization with split tiles and optional
  barplots
- [`gg_seq()`](https://ignophi.github.io/ggrecipes/reference/gg_seq.html) -
  Sequence coverage plot with region highlighting
- [`gg_seqdiff()`](https://ignophi.github.io/ggrecipes/reference/gg_seqdiff.html) -
  Sequence alignment showing only differences
- [`gg_biodist()`](https://ignophi.github.io/ggrecipes/reference/gg_biodist.html) -
  Biodistribution plots with easy free-scale faceting

**Chemoinformatics:**

Binding kinetics and drug discovery visualizations

- [`gg_kdmap()`](https://ignophi.github.io/ggrecipes/reference/gg_kdmap.html) -
  Kinetic rate maps (association/dissociation) with iso-affinity lines

## Installation

Install from CRAN:

``` r
# currently unavailable - initial submission in progress
install.packages("ggrecipes")
```

Development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("Ignophi/ggrecipes")
```

## Documentation

- [Tutorial](https://ignophi.github.io/ggrecipes/articles/ggrecipes.html)
- [Functions](https://ignophi.github.io/ggrecipes/reference/)

## Citation

If you use ggrecipes in your work, please cite:

    [Citation information will be added upon publication]

## License

MIT License - see LICENSE file for details.
