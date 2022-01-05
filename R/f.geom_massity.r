
#' Densité non normalisée
#'
#' Affiche une fonction de densité, discrétisée, exprimée en unité de y de sorte que la surface est égale à la somme de y sur l'espace des x
#' procure quelques fonction accessibles par after_stat (density, mass, cummass, totalmass, mean)
#' permet de comparer des "densités" entre elles puisqu'elles sont exprimées en unité commune (et non pas normées à 1 en surface)
#'
#' @seealso
#'   [geom_massogram()] densité en masse, en transformant l'échelle des x pour avoir une distribution par quantile,
#'   [geom_quantogram()] calcule la densité et normalise les valeurs de x pour que chaque intervalle de x soit un quantile de la distribution de x,
#'   [geom_quantagram()] calcule la densité et normalise les valeurs de x mais conserve l'échelle initiale de x
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_bar
#' @param outline.type encadrement des rectangles ("both", "upper", "lower", "full")
#'
#' @return une "layer" pour ggplot
#' @export
#' @import ggplot2
#'
geom_massity <- function(mapping = NULL, data = NULL,
                         stat = "massity", position = "identity",
                         ...,
                         na.rm = FALSE,
                         orientation = NA,
                         show.legend = NA,
                         inherit.aes = TRUE,
                         outline.type = "upper") {
  outline.type <- match.arg(outline.type, c("both", "upper", "lower", "full"))

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomMassity,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      orientation = orientation,
      outline.type = outline.type,
      ...
    )
  )
}

#' @rdname geom_massity
#' @format NULL
#' @usage NULL
#'
#' @export
GeomMassity <- ggplot2::ggproto("GeomMassity", ggplot2::GeomArea,
                       default_aes = plyr::defaults(
                         ggplot2::aes(fill = NA, mass = 1, colour = "black", alpha = NA),
                         ggplot2::GeomArea$default_aes
                       )
)

#' @rdname geom_massity
#' @format NULL
#' @usage NULL
#'
#' @export
stat_massity <- function(mapping = NULL, data = NULL,
                         geom = "area", position = "stack",
                         ...,
                         bw = "SJ",
                         adjust = 1,
                         kernel = "gaussian",
                         n = 512,
                         trim = FALSE,
                         na.rm = FALSE,
                         orientation = NA,
                         show.legend = NA,
                         inherit.aes = TRUE) {

  layer(
    data = data,
    mapping = mapping,
    stat = StatMassity,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bw = bw,
      adjust = adjust,
      kernel = kernel,
      n = n,
      trim = trim,
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  )
}

#' @rdname geom_massity
#' @format NULL
#' @usage NULL
#' @export
StatMassity <- ggplot2::ggproto("StatMassity", ggplot2::Stat,
                       required_aes = "x|y",

                       default_aes = ggplot2::aes(x = ggplot2::after_stat(mass),
                                         y = ggplot2::after_stat(mass),
                                         fill = NA, mass = NULL),

                       setup_params = function(data, params) {
                         params$flipped_aes <- ggplot2::has_flipped_aes(data, params, main_is_orthogonal = FALSE, main_is_continuous = TRUE)

                         has_x <- !(is.null(data$x) && is.null(params$x))
                         has_y <- !(is.null(data$y) && is.null(params$y))
                         if (!has_x && !has_y) {
                           abort("stat_massity() requires an x or y aesthetic.")
                         }

                         params
                       },

                       extra_params = c("na.rm", "orientation"),

                       compute_group = function(data, scales, bw = "nrd0", adjust = 1, kernel = "gaussian",
                                                n = 512, trim = FALSE, na.rm = FALSE, flipped_aes = FALSE) {
                         data <- ggplot2::flip_data(data, flipped_aes)
                         if (trim) {
                           range <- range(data$x, na.rm = TRUE)
                         } else {
                           range <- scales[[flipped_names(flipped_aes)$x]]$dimension()
                         }

                         density <- compute_massity(data$x, data$mass, from = range[1],
                                                    to = range[2], bw = bw, adjust = adjust, kernel = kernel, n = n)
                         density$flipped_aes <- flipped_aes
                         flip_data(density, flipped_aes)
                       }

)

# compute massity -------------------
# fonction qui calcule les densités, utilisée dans StatMassity

compute_massity <- function(x, m, from, to, bw = "nrd0", adjust = 1,
                            kernel = "gaussian", n = 512) {
  nax <- is.na(x)
  x <- x[!nax]
  nx <- length(x)
  if (is.null(m)) {
    w <- rep(1 / nx, nx)
    tm <- 1
  } else {
    nam <- is.na(m)
    x <- x[!nam]
    nx <- length(x)
    m <- m[!nam]
    tm <- sum(m)
    if(tm>0)
      w <- m / tm
    else
      w <- 1
  }

  # if less than 2 points return data frame of NAs and a warning
  if (nx < 2 | tm==0) {
    warn("Groups with fewer than two data points or 0 mass have been dropped.")
    return(vctrs::new_data_frame(list(
      x = NA_real_,
      density = NA_real_,
      mass = NA_real_,
      mean = NA_real_,
      cummass = NA_real_,
      totalmass = NA_real_,
      n = NA_integer_
    ), n = 1L))
  }

  dens <- stats::density(x, weights = w, bw = bw, adjust = adjust,
                         kernel = kernel, n = n, from = from, to = to)
  uw_dens <- stats::density(x, bw = bw, adjust = adjust,
                            kernel = kernel, n = n, from = from, to = to)

  vctrs::new_data_frame(list(
    x = dens$x,
    density = dens$y,
    mass =  dens$y * tm,
    mean = dens$y * tm / (uw_dens$y*nx),
    cummass = cumsum(dens$y*tm)*(max(dens$x)-min(dens$x))/n,
    totalmass = tm,
    n = nx
  ), n = length(dens$x))
}
