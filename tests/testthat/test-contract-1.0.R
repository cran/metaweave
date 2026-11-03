test_that("constructors reject ambiguous species and invalid probabilities", {
  expect_error(new_assemblage(list(a = c("x", "x"))), "unique")
  expect_error(new_assemblage(setNames(list("x", "y"), c("a", "a"))), "unique")
  expect_error(new_assemblage(list(a = NA_character_)), "non-missing")
  expect_error(new_assemblage(list(a = "x"), cell_id = 0), "cell_id")
  expect_error(new_assemblage(list(a = "x"), coordinates = c(1, NA)), "coordinates")
  expect_error(new_ecological_network(matrix(NA_real_, 1), "probability"), "finite")
  expect_error(new_ecological_network(matrix(2, 1), "binary"), "0 or 1")
  expect_error(new_ecological_network(matrix(-1, 1), "weighted"), "non-negative")
})

test_that("block model validates its complete probability lookup", {
  rows <- data.frame(species = "p", guild = "plant")
  cols <- data.frame(species = "a", guild = "animal")
  theta <- matrix(.5, 1, 1, dimnames = list("plant", "animal"))
  expect_s3_class(block_model(rows, cols, theta, "plants", "animals"), "block_model")
  expect_error(block_model(rows, cols, matrix(NA_real_, 1, 1), "plants", "animals"), "theta")
  expect_error(block_model(rows, cols, matrix(.5, 1, 1,
                                              dimnames = list("wrong", "animal")),
                           "plants", "animals"), "Lookup guilds")
})

test_that("empty collection keeps a typed index and summaries", {
  model <- new_inference_model(function(model, assemblage) {
    new_ecological_network(matrix(.5, 1, 1), "probability")
  })
  collection <- infer_networks(list(), model)
  expect_s3_class(collection, "network_collection")
  expect_equal(names(collection$index), c("cell_id", "x", "y"))
  expect_equal(nrow(collection$index), 0L)
  summarized <- summarize_networks(collection)
  expect_equal(names(summarized$index), c("cell_id", "x", "y",
                                           "row_richness", "column_richness",
                                           "possible_links", "expected_links",
                                           "mean_probability"))
  expect_equal(nrow(summarized$index), 0L)
})

test_that("spatial workflow retains shape with no eligible cells", {
  grid <- create_standard_grid(c(0, 2, 0, 2), resolution = 1)
  plants <- terra::rast(grid)
  animals <- terra::rast(grid)
  names(plants) <- "p"
  names(animals) <- "a"
  terra::values(plants) <- 0
  terra::values(animals) <- 0
  model <- probability_matrix_model(matrix(.5, 1, 1,
    dimnames = list("p", "a")), "plants", "animals")
  result <- run_spatial_inference(list(plants = plants, animals = animals), model)
  expect_equal(length(result$result$networks), 0L)
  expect_equal(nrow(result$result$index), 0L)
  expect_s3_class(result$spatial, "spatial_result")
  expect_equal(names(result$spatial$data),
               c("row_richness", "column_richness", "possible_links",
                 "expected_links", "mean_probability"))
  expect_true(all(is.na(terra::values(result$spatial$data))))
})

test_that("spatial input rejects duplicate names and negative values", {
  grid <- create_standard_grid(c(0, 1, 0, 1), resolution = 1)
  names(grid) <- "a"
  terra::values(grid) <- -1
  expect_error(assemble_communities(list(one = grid)), "non-negative")
  terra::values(grid) <- 1
  expect_error(assemble_communities(setNames(list(grid, grid), c("one", "one"))),
               "unique group")
  expect_error(assemble_communities(list(one = grid), min_species = 0),
               "positive integer")
})

test_that("model outputs report local species excluded by matching", {
  p <- matrix(.7, 1, 1, dimnames = list("plant_a", "animal_a"))
  site <- new_assemblage(list(plants = c("plant_a", "plant_extra"),
                             animals = c("animal_a", "animal_extra")))
  pair <- infer_network(site, probability_matrix_model(p, "plants", "animals"))
  expect_identical(pair$metadata$excluded_rows, "plant_extra")
  expect_identical(pair$metadata$excluded_columns, "animal_extra")

  rows <- data.frame(species = "plant_a", guild = "plant_guild")
  cols <- data.frame(species = "animal_a", guild = "animal_guild")
  theta <- matrix(.7, 1, 1, dimnames = list("plant_guild", "animal_guild"))
  block <- infer_network(site, block_model(rows, cols, theta, "plants", "animals"))
  expect_identical(block$metadata$excluded_rows, "plant_extra")
  expect_identical(block$metadata$excluded_columns, "animal_extra")
})

test_that("count simulation validates inputs and preserves seeded RNG state", {
  p <- matrix(c(.25, .75), 1, 2)
  set.seed(91)
  before <- .Random.seed
  one <- simulate_interaction_counts(p, size = 20, seed = 42)
  expect_identical(.Random.seed, before)
  expect_equal(one, simulate_interaction_counts(p, size = 20, seed = 42))
  expect_equal(sum(one), 20)
  expect_error(simulate_interaction_counts(p, size = -1), "size")
  expect_error(simulate_interaction_counts(p, size = 1.5), "size")
  expect_error(simulate_interaction_counts(matrix(-.1, 1), size = 2), "nonnegative")
  expect_error(simulate_interaction_counts(p, seed = NA_real_), "seed")
  expect_error(simulate_network(p), "ecological network")
})
