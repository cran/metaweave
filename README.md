# metaweave

Spatial ecological network inference from species distributions and supplied
interaction models. Supports species-pair probabilities, stochastic block models,
and directed constrained network ensembles.

`maxent_model()` is experimental. Its fixed-link route samples uniformly under
an exact link-count constraint; its degree-sequence route uses approximate
swap-chain sampling without a built-in convergence diagnostic. It is outside
the stable 1.x API commitment.

The 1.0 API keeps the core spatial workflow, supported specialist helpers, and
existing compatibility wrappers stable through 1.x. `rarefied_network_metric()`
is also experimental while its metric behavior is evaluated.

## Installation

Install the published release from CRAN:

```r
install.packages("metaweave")
```

## A minimal example

```r
library(metaweave)
p <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE,
            dimnames = list(c("plant_a", "plant_b"), c("animal_x", "animal_y")))
model <- probability_matrix_model(p, "plants", "animals")
site <- new_assemblage(list(plants = c("plant_a", "plant_b"), animals = "animal_x"))
network <- infer_network(site, model)
summarize_network(network)
network$metadata$excluded_rows
```

For a complete runnable spatial example, use
`vignette("getting-started", package = "metaweave")`.
The tutorial uses synthetic data; publication case studies are distributed
separately. Use `citation("metaweave")` for the software citation.

Species present locally but absent from a supplied pair matrix or block lookup
are excluded from the inferred network. Inspect `network$metadata$excluded_rows`
and `network$metadata$excluded_columns` to find them. MetaWeave uses supplied
interaction probabilities; it does not estimate or validate their ecological
calibration.

Author, maintainer and copyright holder: Gabriel Munoz, Concordia University.
Contact: gmunozacevedo@proton.me. License: MIT.
