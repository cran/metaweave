## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(metaweave)
probabilities <- matrix(
  c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE,
  dimnames = list(c("plant_a", "plant_b"), c("animal_x", "animal_y"))
)
model <- probability_matrix_model(probabilities, "plants", "animals")

## -----------------------------------------------------------------------------
site <- new_assemblage(list(plants = c("plant_a", "plant_b"),
                            animals = "animal_x"))
network <- infer_network(site, model)
network$matrix
summarize_network(network)

## -----------------------------------------------------------------------------
grid <- create_standard_grid(c(0, 2, 0, 2), resolution = 1)
plants <- terra::rast(list(grid, grid))
animals <- terra::rast(list(grid, grid))
names(plants) <- rownames(probabilities)
names(animals) <- colnames(probabilities)
terra::values(plants) <- cbind(c(1, 1, 0, 1), c(1, 0, 1, 1))
terra::values(animals) <- cbind(c(1, 1, 1, 0), c(0, 1, 1, 1))
result <- run_spatial_inference(
  distributions = list(plants = plants, animals = animals),
  model = model, min_species = 1
)
result$result$index

## ----map, fig.width=5, fig.height=4-------------------------------------------
terra::plot(result$spatial$data[["expected_links"]], main = "Expected links per cell")

