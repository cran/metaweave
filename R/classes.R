#' Construct a distribution collection
#'
#' @param data A nonempty `terra::SpatRaster` with unique, nonempty species layer names.
#' @param guild Optional guild label for the collection.
#' @param metadata Optional species metadata data frame.
#' @return A `distribution_collection` list with `data` (the raster), `guild`, and `metadata`.
#' @export
new_distribution_collection <- function(data, guild = NULL, metadata = NULL) {
  if (!inherits(data, "SpatRaster")) stop("`data` must be a SpatRaster.")
  if (!length(data) || is.null(names(data)) || anyNA(names(data)) ||
      any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("Raster layers must have unique, non-empty species names.")
  }
  structure(list(data = data, guild = guild, metadata = metadata), class = "distribution_collection")
}

#' @export
print.distribution_collection <- function(x, ...) {
  cat("<distribution_collection>\n", length(x$data), " species\n", sep = "")
  invisible(x)
}

#' Construct a local assemblage
#' @param species Nonempty named list of species vectors, one per ecological
#'   group. Group names and species names within a group must be unique and
#'   nonempty; species names are coerced to character.
#' @param cell_id Optional positive integer raster cell identifier or `NA`.
#' @param coordinates Optional two finite numeric values in x/y order.
#' @return An `assemblage` list with `species`, `cell_id`, and `coordinates`.
#' @examples
#' new_assemblage(list(plants = c("a", "b"), animals = "x"), cell_id = 1)
#' @export
new_assemblage <- function(species, cell_id = NA_integer_, coordinates = NULL) {
  if (!is.list(species) || !length(species) || is.null(names(species)) ||
      anyNA(names(species)) || any(!nzchar(names(species))) ||
      anyDuplicated(names(species))) {
    stop("`species` must be a non-empty list with unique, non-empty group names.")
  }
  species <- lapply(species, as.character)
  if (any(vapply(species, function(x) anyNA(x) || any(!nzchar(x)) ||
                 anyDuplicated(x) > 0L, logical(1)))) {
    stop("Species names within each group must be unique, non-missing, and non-empty.")
  }
  if (length(cell_id) != 1L || (!is.na(cell_id) &&
      (!is.numeric(cell_id) || !is.finite(cell_id) || cell_id < 1 ||
       cell_id != floor(cell_id)))) {
    stop("`cell_id` must be one positive integer or NA.")
  }
  if (!is.null(coordinates) &&
      (!is.numeric(coordinates) || length(coordinates) != 2L ||
       any(!is.finite(coordinates)))) {
    stop("`coordinates` must contain two finite numeric values (x, y).")
  }
  structure(list(species = species, cell_id = cell_id, coordinates = coordinates), class = "assemblage")
}

#' @export
print.assemblage <- function(x, ...) {
  cat("<assemblage>", if (!is.na(x$cell_id)) paste(" cell", x$cell_id), "\n")
  for (nm in names(x$species)) cat(nm, ": ", length(x$species[[nm]]), " species\n", sep = "")
  invisible(x)
}

#' Construct a pluggable inference model
#' @param predict A function accepting `(model, assemblage)` and returning an ecological network.
#' @param name Model name.
#' @param parameters Model-specific parameters.
#' @param subclass Optional additional class name.
#' @return An `inference_model` list containing `name`, `predict`, and `parameters`, with the optional subclass prepended.
#' @export
new_inference_model <- function(predict, name = "custom", parameters = list(), subclass = NULL) {
  if (!is.function(predict)) stop("`predict` must be a function.")
  structure(list(name = name, predict = predict, parameters = parameters),
            class = c(subclass, "inference_model"))
}

#' @export
print.inference_model <- function(x, ...) {
  cat("<inference_model>", x$name, "\n")
  invisible(x)
}

#' Construct an ecological network
#' @param matrix Numeric interaction matrix. Probabilities must be finite in
#'   `[0, 1]`; binary entries must be 0 or 1. Weighted and observed entries
#'   must be nonnegative and may be `NA`.
#' @param type One of probability, binary, weighted, or observed.
#' @param metadata Optional metadata.
#' @param row_group,column_group Ecological groups represented by rows and columns.
#' @param directed Whether matrix orientation represents directed interactions.
#' @param self_links Whether self-links are permitted; `NA` when not applicable.
#' @return An `ecological_network` list with `matrix`, `type`, `metadata`, `row_group`, `column_group`, `directed`, and `self_links`; also inherits from the type-specific network class.
#' @export
new_ecological_network <- function(matrix, type = c("probability", "binary", "weighted", "observed"), metadata = list(),
                                   row_group = "rows", column_group = "columns",
                                   directed = FALSE, self_links = NA) {
  type <- match.arg(type)
  if (!is.matrix(matrix) || !is.numeric(matrix)) stop("`matrix` must be a numeric matrix.")
  if (type == "probability" &&
      (any(!is.finite(matrix)) || any(matrix < 0 | matrix > 1))) {
    stop("Probabilities must be finite and lie in [0, 1].")
  }
  if (type == "binary" && (any(!is.finite(matrix)) ||
      any(!(matrix %in% c(0, 1))))) {
    stop("Binary network entries must be 0 or 1.")
  }
  if (type %in% c("weighted", "observed") &&
      (any(is.infinite(matrix) | is.nan(matrix)) ||
       any(matrix < 0, na.rm = TRUE))) {
    stop("Network weights must be non-negative and not infinite.")
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) stop("`directed` must be TRUE or FALSE.")
  if (!is.logical(self_links) || length(self_links) != 1L) stop("`self_links` must be TRUE, FALSE, or NA.")
  structure(list(matrix = matrix, type = type, metadata = metadata,
                 row_group = as.character(row_group), column_group = as.character(column_group),
                 directed = directed, self_links = self_links),
            class = c(paste0(type, "_network"), "ecological_network"))
}

#' @export
print.ecological_network <- function(x, ...) {
  cat("<ecological_network:", x$type, ">", nrow(x$matrix), x$row_group,
      "x", ncol(x$matrix), x$column_group,
      if (isTRUE(x$directed)) "directed" else "undirected", "\n")
  invisible(x)
}

#' Construct a collection of spatial ecological networks
#' @param networks List of ecological networks.
#' @param index Cell-level data frame.
#' @param template Spatial raster template.
#' @param settings Optional workflow settings.
#' @return A `network_collection` list with `networks`, `index`, `template`, and `settings`. Each index row corresponds to one network.
#' @export
new_network_collection <- function(networks, index, template = NULL, settings = list()) {
  if (!is.list(networks) || !is.data.frame(index) || length(networks) != nrow(index)) {
    stop("`networks` and rows of `index` must have equal length.")
  }
  structure(list(networks = networks, index = index, template = template, settings = settings),
            class = "network_collection")
}

#' @export
print.network_collection <- function(x, ...) {
  cat("<network_collection>", length(x$networks), " local networks\n")
  invisible(x)
}

#' Construct a spatial result
#' @param data A spatial object.
#' @param metrics Names of mapped metrics.
#' @return A `spatial_result` list with the spatial object in `data` and metric names in `metrics`.
#' @export
new_spatial_result <- function(data, metrics = names(data)) {
  structure(list(data = data, metrics = metrics), class = "spatial_result")
}
