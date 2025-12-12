# ggrecipes <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ggrecipes)](https://CRAN.R-project.org/package=ggrecipes)
<!-- badges: end -->

A collection of easy-to-use custom [ggplot2](https://ggplot2.tidyverse.org/index.html)-based functions for data exploration and analysis. Each function handles data preprocessing and returns a ggplot2 object that can be further customized using standard ggplot2 syntax. Includes general-purpose and domain-specific visualizations.

## Features

![](man/figures/demo.gif)

**General:**

General-purpose visualizations

- `gg_splitcorr()` - Split-Correlation Heatmap 
- `gg_rankshift()` - Paired box-/bar- plots with rank change
- `gg_criteria()` - Criteria heatmap with optional barplots
- `gg_conf()` - Confusion/contingency table bubble plot

**Bioinformatics:**

Sequence analysis and genomics visualizations

- `gg_geno()` - Biallelic genotype visualization with split tiles and optional barplots
- `gg_seq()` - Sequence coverage plot with region highlighting
- `gg_seqdiff()` - Sequence alignment showing only differences
- `gg_biodist()` - Biodistribution plots with easy free-scale faceting

**Chemoinformatics:**

Binding kinetics and drug discovery visualizations

- `gg_kdmap()` - Kinetic rate maps (association/dissociation) with iso-affinity lines


## Installation

Install from CRAN:
```r
# currently unavailable - initial submission in progress
install.packages("ggrecipes")
```

Development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("Ignophi/ggrecipes")
```

## Documentation

- [Function reference](https://ignophi.github.io/ggrecipes/reference/)
- [Vignettes and examples](https://ignophi.github.io/ggrecipes/articles/)

## Citation

If you use ggrecipes in your work, please cite:
```
[Citation information will be added upon publication]
```

## License

MIT License - see LICENSE file for details.
