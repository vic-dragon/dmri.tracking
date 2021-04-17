
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dmri.tracking

<!-- badges: start -->

<!-- badges: end -->

The goal of dmri.tracking is to apply the deterministic tracking
algorithm from DiST (Raymond et al. 2016) and to visualize the tracking
result.

## Installation

You can install the released version of dmri.tracking from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dmri.tracking")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vic-dragon/dmri.tracking")
```

## Example

This is an example of the usage of this library:

    library(dmri.tracking)
    #Load example output from peak detection algorithm
    load(system.file("extdata", "peakresult.rda", package = "dmri.tracking"))
    
    peak.result  #Output from the peak detection algorithm
    
    #Apply Tracking algorithm
    result = v.track(v.obj = peak.result, max.line=500)
    
    library(rgl)
    open3d()
    for (iind in (result$sorted.iinds[result$sorted.update.ind])){
      cat(iind,"\n")
      tractography(result$tracks1[[iind]]$inloc, result$tracks1[[iind]]$dir)
      tractography(result$tracks2[[iind]]$inloc, result$tracks2[[iind]]$dir)
    }
