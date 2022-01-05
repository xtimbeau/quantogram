#' valeur moyenne pondérée avec transformation de l'échelle des x
#'
#' Affiche la valeur moyenne (et les quantiles éventuellement) d'une variable
#' L'échelle des x est transformée en quantiles (et les moyennes sont calculées sur ces bins).
#' procure quelques valeurs accessibles par after_stat (sum_gross, mass)
#' permet de comparer des "valeurs" entre elles sur des espaces normalisés
#'
#' @seealso
#'   [geom_massity()] densité en masse,
#'   [geom_massogram()] ,
#'   [geom_quantagram()] calcule la moyenne pondérée sur des intervalles de x normalisé mais conserve l'échelle initiale de x
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_bar
#' @param outline.type encadrement des rectangles ("both", "upper", "lower", "full")
#'
#' @return une "layer" pour ggplot
#' @export
#'
#' @import ggplot2 dplyr rlang data.table vctrs

geom_quantogram <- function(mapping = NULL, data = NULL,
                            stat = StatQuanto, position = "identity",
                            ...,
                            na.rm = FALSE,
                            cuts = NULL,
                            probs = 0.5,
                            lines = TRUE,
                            bars = TRUE,
                            bins = 0,
                            trans=FALSE,
                            se = TRUE,
                            show.legend = NA,
                            inherit.aes = TRUE)
{
  cuts_trans <- NULL
  if (!is.null(mapping$x) & is_function(cuts)) {
    xx <- rlang::quo_name(mapping$x)
    fcuts <- cuts(xx)
    if (!all(is.na(fcuts))) {
      cuts <- fcuts
      ncuts <- names(fcuts)
      icuts <- stringr::str_extract(ncuts, pattern = "[:digit:]*\\.?[:digit:]*") |> as.numeric()
      cuts_trans <- scales::trans_new(
        "qtrans",
        approxfun(y=icuts, x=fcuts),
        approxfun(x=icuts, y=fcuts),
        breaks = scales::extended_breaks(),
        minor_breaks = scales::regular_minor_breaks(),
        format = scales::format_format(),
        domain = c(min(fcuts), max(fcuts)))
    } else {
      cuts <- NULL
    }
  }
  params <- list(
    na.rm = na.rm,
    bins = bins,
    cuts = cuts,
    probs = sort(probs),
    lines = lines,
    bars = bars,
    trans= trans,
    se=se,
    ...
  )

  if(!is.null(cuts)) {

  }

  list(
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomQuantogram,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = params
    ),
    if(!is.null(cuts_trans)&trans) coord_trans(x=cuts_trans))
}

#' @rdname geom_quantogram
#' @format NULL
#' @usage NULL
#'
#' @export
GeomQuantogram <- ggplot2::ggproto(
  "GeomQuantogram",
  ggplot2::Geom,
  setup_params = function(data, params) {
    params
  },
  required_aes = c("x", "y"),
  optional_aes = c("label",  "w"),
  non_missing_aes = c("ymin", "ymax", "xend", "yend"),
  default_aes = plyr::defaults(
    ggplot2::aes(fill = NA, colour = "black", alpha = NA, size = 0.5),
    ggplot2::GeomArea$default_aes),
  extra_params = c("na.rm", "cuts", "probs", "lines", "bins", "trans",
                   "se", "delta_x", "labels_x"),

  draw_group = function(data, panel_params, coord, probs = 0.5, lines = TRUE, bars = TRUE, se = TRUE, flipped_aes = FALSE) {
    np <- length(probs)
    alphas <- if(np>1)
      c(rep(0.15/(np-1), np-1), 0.1)
    else
      0.25

    if(se){
      if (lines) {
        ribbons <- purrr::map2(probs, alphas, ~ data |>
                          dplyr::transmute(x,
                                    PANEL,
                                    group,
                                    y = median,
                                    ymin = .data[[glue::glue("y_{.x}_m")]],
                                    ymax = .data[[glue::glue("y_{.x}_p")]],
                                    fill = colour,
                                    colour = NA,
                                    linetype,
                                    size,
                                    alpha = .y
                          ))
      } else {
        ribbons <- purrr::map2(probs, alphas, ~ data |>
                          dplyr::transmute(x,
                                    PANEL,
                                    group,
                                    y = median,
                                    ymin = .data[[glue::glue("y_{.x}_m")]],
                                    ymax = .data[[glue::glue("y_{.x}_p")]],
                                    colour = colour,
                                    linetype,
                                    size,
                                    alpha = .y
                          ))
      }
    }
    else
      ribbons <- NULL

    # line <- data |> transmute(colour, x, y, ymin = 0, ymax = y, PANEL, group, fill, size, linetype, alpha)

    line <- data |>
      dplyr::transmute(colour, x=x-dx/2, xend=x+dx, y, yend=y, PANEL, group, fill, size, linetype, alpha)
    line2 <- data |>
      dplyr::arrange(x) |>
      dplyr::transmute(colour, x=x+dx/2, xend=x, y, yend=lead(y), PANEL, group, fill, size, linetype, alpha)
    line <- bind_rows(line, line2) |> drop_na(x,xend,y,yend)

    rects <- data |>
      dplyr::arrange(x) |>
      dplyr::group_by(PANEL, group) |>
      dplyr::transmute(colour=NA, xmin=if_else(is.na(lag(x)), x, x-(x-lag(x))/2), xmax=if_else(is.na(lead(x)), x, x+(lead(x)-x)/2),
                ymin = 0, ymax = y, PANEL, group, fill, size, linetype, alpha) |>
      tidyr::drop_na(xmin, xmax) |>
      dplyr::ungroup()
# browser()
    if (lines) {
      ll <- append(
        if (se) purrr::map(ribbons, ~ ggplot2::GeomArea$draw_panel(.x, panel_params, coord)),
        list(ggplot2::GeomSegment$draw_panel(line, panel_params, coord))
      )
    } else {
      ll <- append(
        if (se) purrr::map(ribbons, ~ ggplot2::GeomLinerange$draw_panel(.x, panel_params, coord)),
        list(ggplot2::GeomSegment$draw_panel(line, panel_params, coord),
             ggplot2::GeomRect$draw_panel(rects, panel_params, coord))
      )
    }
    do.call(grid::gList, ll)
  },

  draw_key = ggplot2::draw_key_polygon
)

# stat_quantogram <- function(mapping = NULL, data = NULL,
#                             geom = "area", position = "identity",
#                             ...,
#                             na.rm = FALSE,
#                             cuts = NULL,
#                             probs = c(0.25, 0.75),
#                             lines = TRUE,
#                             bars = TRUE,
#                             trans = FALSE,
#                             show.legend = NA,
#                             inherit.aes = TRUE) {
#   ggplot2::layer(
#     data = data,
#     mapping = mapping,
#     stat = StatQuanto,
#     position=position,
#     geom = geom,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(
#       na.rm = na.rm,
#       bins = bins,
#       cuts = cuts,
#       probs = sort(probs),
#       lines = lines,
#       bars = bars,
#       trans=trans,
#       ...
#     )
#   )
# }

# Fonction pour le compute group: bin les x et calcule les quantiles sur le x biné
#' @rdname geom_quantogram
#' @format NULL
#' @usage NULL
#'
#' @export
StatQuanto <- ggplot2::ggproto(
  "StatQuanto",
  ggplot2::Stat,
  required_aes = c("x", "mass"),
  default_aes = ggplot2::aes(y = ggplot2::after_stat(median),
                    fill = NA, w = 1, alpha = NA, colour = NA, size = 0.5),
  optional_aes = c("label",  "w"),
  non_missing_aes = c("y", "ymin", "ymax"),
  setup_params = function(self, data, params) {
    has_x <- !(is.null(data$x) && is.null(params$x))
    has_measure <- !(is.null(data$mass) && is.null(params$mass))

    if (!has_x && !has_measure) {
      rlang::abort("stat_quanto() requires an x and mass aesthetic.")
    }
    if (is.null(params$label)) params$label <- NULL
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
  extra_params = c("na.rm", "cuts", "probs", "lines", "bins", "trans", "labels_x", "delta_x", "se"),
  setup_data = function(data, params) {
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
        delta_x <- (tail(cuts, -1) - head(cuts, -1))
        delta_x <- c(delta_x[[1]], delta_x)
        data <- data |>
          dplyr::mutate( dx = delta_x[which(x==cuts)])
        return(data)
      }
    }
    else {
      cuts <- params$cuts
    }
    labels_x <- (tail(cuts, -1) + head(cuts, -1)) / 2
    delta_x <- tail(cuts, -1) - head(cuts, -1)
    x_cutted <- findInterval(data$x, cuts, all.inside = TRUE)
    data <- data |>
      dplyr::mutate(
        x_cutted = x_cutted,
        x = labels_x[x_cutted],
        dx = delta_x[x_cutted])

    return(data)
  },
  compute_group = function(data, scales, na.rm = FALSE, bars, probs, trans, labels_x, delta_x) {
    if(!bars) probs <- NULL
    compute_quansity_dt(data$x, data$mass, data$dx, prob = probs,
                        trans = trans, labels_x=labels_x, delta_x=delta_x)
  }
)

#'
#' @import data.table
#'
#'
compute_quansity_dt <- function(x, y, dx, prob, trans, labels_x, delta_x) {
  require(data.table, quietly=TRUE)
  # if(trans)
  #   dx <- rep(sum(dx)/length(dx), length(dx))
  # data <- data.table(x = x, y = y, dx = dx)
  data <- data.table(x = x, y = y, dx = dx)
  data <- merge(data, data.table(x=labels_x, dx=delta_x),
                by=c("x", "dx"),
                all.x=TRUE, all.y=TRUE, sort=TRUE)
  data[is.na(y), y:=0]
  if(trans)
    ydx <- mean(data$dx)
  else
    ydx <- data$dx
  data[, ydx := (ydx)]

  quansity1 <- data[,.(median= median(y),
                       mean = mean(y),
                       n = .N,
                       dx=first(dx)), by=x]
  quansity1[, `:=`(mass = mean)]
  quansities <- map(prob, ~{
    qm <- stringr::str_c("y_", .x, "_m")
    qp <- stringr::str_c("y_", .x, "_p")
    qq <- data[, .(
      q_m = quantile(x = y, probs = c(0.5 - .x / 2)),
      q_p = quantile(x = y, probs = c(0.5 + .x / 2))
    ), by=x]
    setnames(qq, c("q_p", "q_m"), c(qp, qm))
    qq[, x:=NULL]})
  quansities <- do.call(cbind, quansities)
  if(length(prob)>0) quansities[, ':='(ymax=do.call(pmax,.SD), ymin=do.call(pmin,.SD))]
  tibble::as_tibble(cbind(quansity1, quansities))
}

# weighted_quantile_x --------------------------

weighted_quantile_x <- function(x, w = NULL, probs = seq(0, 1, 0.25), na.rm = TRUE,
                                names = TRUE, ...)
{
  if (is.null(w) || length(w)==1) {
    return(quantile(x = x, probs = probs, na.rm = na.rm, names = names, ...))
  }
  stopifnot(length(x) == length(w), all(probs >= 0 & probs <= 1))
  if (isTRUE(na.rm)) {
    ok <- !is.na(x)&!is.na(w)
    x <- x[ok]
    w <- w[ok]
  }
  if (length(x)==0)
    out <- rep(NA, length(probs))
  else
  {
    if (all(w == 0))
      out=rep(0, length(probs))
    else if (length(x) == 1L) {
      out <- rep(x, length(probs))
    }
    else {
      ord <- order(x)
      x <- x[ord]
      w <- w[ord]
      w <- w/sum(w)
      cs <- cumsum(w[-length(w)])
      out <- stepfun(cs, x)(probs)
    }
  }
  if (names) {
    names(out) <- paste0(format(100 * probs, trim = TRUE),"%")
  }
  out
}
