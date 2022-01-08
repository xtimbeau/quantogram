#' Densité non normalisée avec transformation de l'échelle des x
#'
#' Affiche une fonction de densité, discrétisée, exprimée en unité de y de sorte que la surface est égale à la somme de y sur l'espace des x.
#' L'échelle des x est transformée en quantiles pour garantir que la surface de chaque rectangle est bien égale à la somme des y sur cette intervalle.
#' procure quelques fonction accessibles par after_stat (sum_gross, mass)
#' permet de comparer des "densités" entre elles puisqu'elles sont exprimées en unité commune (et non pas normées à 1 en surface)
#'
#' @seealso
#'   [geom_massity()] densité en masse,
#'   [geom_quantogram()] calcule la densité et normalise les valeurs de x pour que chaque intervalle de x soit un quantile de la distribution de x,
#'   [geom_quantagram()] calcule la densité et normalise les valeurs de x mais conserve l'échelle initiale de x
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_bar
#' @param outline.type encadrement des rectangles ("both", "upper", "lower", "full")
#'
#' @return une "layer" pour ggplot
#' @export
#'
#' @import ggplot2
#'

geom_massogram <- function(mapping = NULL, data = NULL,
                            stat = StatMasso,
                            position = "stack",
                            ...,
                            na.rm = FALSE,
                            cuts = NULL,
                            lines = FALSE,
                            bins = 0,
                            trans=FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE)
{
  cuts_trans <- NULL
  if (!is.null(mapping$x) & is_function(cuts)) {
    xx <- rlang::quo_name(mapping$x)
    fcuts <- cuts(xx)
    if (!all(is.na(fcuts))) {
      cuts <- fcuts
      icuts <- stringr::str_extract(names(fcuts), pattern = "[:digit:]*\\.?[:digit:]*") |> as.numeric()
      cuts_trans <- scales::trans_new(
        "qtrans",
        approxfun(y=icuts, x=fcuts, rule=2),
        approxfun(x=icuts, y=fcuts, rule=2),
        breaks = scales::extended_breaks(),
        minor_breaks = scales::regular_minor_breaks(),
        format = scales::format_format(),
        domain = c(min(fcuts), max(fcuts)))
    }
  }
  params <- list(
    na.rm = na.rm,
    bins = bins,
    cuts = cuts,
    lines = lines,
    trans= trans,
    ...
  )

  if(!is.null(cuts)) {

  }

  list(
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomMassogram,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = params
    ),
    if(!is.null(cuts_trans)&trans) ggplot2::coord_trans(x=cuts_trans))
}

#' @rdname geom_massogram
#' @format NULL
#' @usage NULL
#'
#' @export
GeomMassogram <- ggplot2::ggproto(
  "GeomMassogram",
  ggplot2::GeomArea,
  setup_params = function(data, params) {
    params
  },

  extra_params = c("na.rm", "cuts", "lines", "bins", "trans"),

  draw_group = function(self, data, panel_params, coord, lines) {
    if (lines)
    {
      line <- data |>
        dplyr::transmute(colour, x, y, ymin=ymin, ymax=y, PANEL, group, fill, size, linetype, alpha)
      ll <- ggplot2::GeomArea$draw_panel(line, panel_params, coord)
    }
    else
    {
      line <- data |>
        dplyr::transmute(colour, x=x-dx/2, xend=x+dx, y, yend=y, PANEL, group, fill, size, linetype, alpha)
      line2 <- data |>
        dplyr::arrange(x) |>
        dplyr::transmute(colour, x=x+dx/2, xend=x, y, yend=lead(y), PANEL, group, fill, size, linetype, alpha)
      line <- dplyr::bind_rows(line, line2)

      rects <- data |>
        dplyr::transmute(colour=NA, xmin=x-dx/2, xmax=x+dx/2,
                  ymin = ymin, ymax = y, PANEL, group, fill, size, linetype, alpha)
      ll <- grid::gList(
        ggplot2::ggproto_parent(GeomRect, self)$draw_panel(rects, panel_params, coord),
        ggplot2::ggproto_parent(GeomSegment, self)$draw_panel(line, panel_params, coord))
    }
    ll
  },

  draw_key = ggplot2::draw_key_polygon,
  required_aes = c("x", "mass"),
  optional_aes = c("label", "w", "xend", "yend", "xmin", "xmax", "ymin", "ymax"),
  default_aes = plyr::defaults(
    ggplot2::aes( fill = NA, colour = "black", alpha = NA, size = 0.5),
    ggplot2::GeomArea$default_aes)
)

#' @rdname geom_massogram
#' @format NULL
#' @usage NULL
#'
#' @export
stat_massogram <- function(mapping = NULL, data = NULL,
                            geom = "area",
                            position = "identity",
                            ...,
                            na.rm = FALSE,
                            cuts = NULL,
                            lines = FALSE,
                            trans = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatMasso,
    position=position,
    geom = geom,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      bins = bins,
      cuts = cuts,
      lines = lines,
      trans=trans,
      ...
    )
  )
}

# Fonction pour le compute group: bin les x et calcule les quantiles sur le x biné -----------
#' @rdname geom_massogram
#' @format NULL
#' @usage NULL
#'
#' @export
StatMasso <- ggplot2::ggproto(
  "StatMasso",
  ggplot2::Stat,
  required_aes = c("x", "mass"),
  default_aes = ggplot2::aes(y = ggplot2::after_stat(sum), fill = NA, w = 1, alpha = NA, colour = NA, size = 0.5),
  optional_aes = c("label", "ymin", "ymax", "w"),

  setup_params = function(self, data, params) {
    has_x <- !(is.null(data$x) && is.null(params$x))
    has_mass <- !(is.null(data$mass) && is.null(params$mass))

    if (!has_x && !has_mass) {
      abort("stat_massogram() requires an x and mass aesthetic.")
    }
    has_cuts <- !is.null(params$cuts)
    has_weights <- !is.null(data$w)
    has_label <- !is.null(data$label)

    if (!has_cuts) {
      if (has_label) {
        xx <- tibble::tibble(x = data$x, label = data$label) |>
          dplyr::group_by(label) |>
          dplyr::summarize(x = median(x))
      }
      else {
        xx <- tibble::tibble(x = data$x)
      }
      if (has_weights) {
        w <- data$w
      } else {
        w <- 1
      }
      if (params$bins > 0) {
        cuts <- weighted_quantile_x(xx$x, w = w, prob = 0:params$bins / params$bins)
      } else {
        cuts <- sort(unique(xx$x))
      }
    }
    else {
      cuts <- params$cuts
    }
    delta_x <- (tail(cuts, -1) - head(cuts, -1))
    params$cuts <- cuts
    params$delta_x <- delta_x
    params$labels_x <- (tail(cuts, -1) + head(cuts, -1)) / 2
    params
  },
  extra_params = c("na.rm", "cuts", "lines", "bins", "trans", "labels_x", "delta_x"),
  setup_data = function(data, params) {
    x_cutted <- findInterval(data$x, params$cuts, all.inside = TRUE)
    data <- data |>
      dplyr::mutate(
        x_cutted = x_cutted,
        x = params$labels_x[x_cutted],
        dx = params$delta_x[x_cutted],
        gm = sum(mass)) |>
      dplyr::group_by(PANEL, group) |>
      dplyr::mutate(ggm=sum(mass, na.rm=TRUE)) |>
      dplyr::ungroup()
    return(data)
  },
  compute_group = function(data, scales, na.rm = FALSE, probs, trans, labels_x, delta_x) {
    compute_masso(data$x, data$mass, data$dx, data$ggm, trans = trans, labels_x=labels_x, delta_x=delta_x)
  }
)

# fonction interne pour le calcul

#'
#' @import data.table

compute_masso <- function(x, y, dx, ggm, trans=FALSE, labels_x, delta_x) {
  data <- data.table(x = x, y = y, dx = dx)
  data <- merge(data,
                data.table(x=labels_x, dx=delta_x),
                by=c("x", "dx"), all.x=TRUE, all.y=TRUE, sort=TRUE)
  data[is.na(y), `:=`(y=0)]

  if(trans)
    ydx <- mean(data$dx)
  else
    ydx <- data$dx
  data[, `:=`(ydx = (ydx))]

  quansity1 <- data[,.(sum_gross = sum(y),
                       mass = sum(y),
                       dx=first(dx),
                       ydx=first(ydx)), by=x]
  setorder(quansity1, x)
  quansity1[, `:=`(sum = sum_gross/ydx)][, `:=`(cumsum = cumsum(sum_gross))]
  quansity1[, `:=`(ymax = sum, ymin=0, groupmass = sum(sum_gross), grossgroupmass = ggm[[1]])]
  tibble::as_tibble(quansity1)

}
