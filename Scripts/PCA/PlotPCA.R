library(ggfortify)

my_plot_label = function (p, data, x = NULL, y = NULL, label = TRUE, label.label = "rownames", 
                          label.colour = NULL, label.alpha = NULL, label.size = NULL, 
                          label.angle = NULL, label.family = NULL, label.fontface = NULL, 
                          label.lineheight = NULL, label.hjust = NULL, label.vjust = NULL, 
                          label.repel = FALSE, label.show.legend = NA, label.parse = TRUE) 
{
  if (!is.data.frame(data)) {
    stop(paste0("Unsupported class: ", class(data)))
  }
  if (!missing(label.colour) && !is.null(label.colour) && 
      missing(label)) {
    label <- TRUE
  }
  if (label || label.repel) {
    if (is.null(label.colour)) {
      label.colour <- "#000000"
    }
    if (label.repel && "ggrepel" %in% rownames(installed.packages())) {
      textfunc <- ggrepel::geom_text_repel
    }
    else {
      textfunc <- ggplot2::geom_text
    }
    p <- p + ggfortify:::geom_factory(textfunc, data, x = x, y = y, 
                                      label = label.label, colour = label.colour, alpha = label.alpha, 
                                      size = label.size, angle = label.angle, family = label.family, 
                                      fontface = label.fontface, lineheight = label.lineheight, 
                                      hjust = label.hjust, vjust = label.vjust, parse = label.parse,
                                      show.legend = label.show.legend)
  }
  p
}


my_ggbiplot = function (plot.data, loadings.data = NULL, colour = NULL, size = NULL, 
                        linetype = NULL, alpha = NULL, fill = NULL, shape = NULL, 
                        label = FALSE, label.label = "rownames", label.colour = colour, 
                        label.alpha = NULL, label.size = NULL, label.angle = NULL, 
                        label.family = NULL, label.fontface = NULL, label.lineheight = NULL, 
                        label.hjust = NULL, label.vjust = NULL, label.repel = FALSE, 
                        label.parse = TRUE,
                        loadings = FALSE, loadings.colour = "#FF0000", loadings.label = FALSE, 
                        loadings.label.label = "rownames", loadings.label.colour = "#FF0000", 
                        loadings.label.alpha = NULL, loadings.label.size = NULL, 
                        loadings.label.angle = NULL, loadings.label.family = NULL, 
                        loadings.label.fontface = NULL, loadings.label.lineheight = NULL, 
                        loadings.label.hjust = NULL, loadings.label.vjust = NULL, 
                        loadings.label.repel = FALSE, label.show.legend = NA, 
                        loadings.label.parse = TRUE,
                        frame = FALSE, 
                        frame.type = NULL, frame.colour = colour, frame.level = 0.95, 
                        frame.alpha = 0.2, xlim = c(NA, NA), ylim = c(NA, NA), log = "", 
                        main = NULL, xlab = NULL, ylab = NULL, asp = NULL, ...) 
{
  plot.columns <- colnames(plot.data)
  mapping <- ggplot2::aes_string(x = plot.columns[1L], y = plot.columns[2L])
  if (is.logical(shape) && !shape && missing(label)) {
    label <- TRUE
  }
  p <- ggplot2::ggplot(data = plot.data, mapping = mapping)
  if (!is.logical(shape) || shape) {
    p <- p + ggfortify:::geom_factory(ggplot2::geom_point, plot.data, 
                                      colour = colour, size = size, linetype = linetype, 
                                      alpha = alpha, fill = fill, shape = shape)
  }
  if (loadings.label && !loadings) {
    loadings <- TRUE
  }
  if (loadings && !is.null(loadings.data)) {
    scaler <- min(max(abs(plot.data[, 1L]))/max(abs(loadings.data[, 
                                                                  1L])), max(abs(plot.data[, 2L]))/max(abs(loadings.data[, 
                                                                                                                         2L])))
    loadings.columns <- colnames(loadings.data)
    loadings.mapping <- ggplot2::aes_string(x = 0, y = 0, 
                                            xend = loadings.columns[1L], yend = loadings.columns[2L])
    loadings.data[, 1L:2L] <- loadings.data[, 1L:2L] * scaler * 
      0.8
    p <- p + geom_segment(data = loadings.data, mapping = loadings.mapping, 
                          arrow = grid::arrow(length = grid::unit(8, "points")), 
                          colour = loadings.colour)
    p <- my_plot_label(p = p, data = loadings.data, label = loadings.label, 
                       label.label = loadings.label.label, label.colour = loadings.label.colour, 
                       label.alpha = loadings.label.alpha, label.size = loadings.label.size, 
                       label.angle = loadings.label.angle, label.family = loadings.label.family, 
                       label.fontface = loadings.label.fontface, label.lineheight = loadings.label.lineheight, 
                       label.hjust = loadings.label.hjust, label.vjust = loadings.label.vjust, 
                       label.repel = loadings.label.repel, label.show.legend = label.show.legend,
                       label.parse = loadings.label.parse)
  }
  p <- my_plot_label(p = p, data = plot.data, label = label, 
                     label.label = label.label, label.colour = label.colour, 
                     label.alpha = label.alpha, label.size = label.size, 
                     label.angle = label.angle, label.family = label.family, 
                     label.fontface = label.fontface, label.lineheight = label.lineheight, 
                     label.hjust = label.hjust, label.vjust = label.vjust, 
                     label.repel = label.repel, label.show.legend = label.show.legend,
                     label.parse = label.parse)
  if (missing(frame) && !is.null(frame.type)) {
    frame <- TRUE
  }
  . <- NULL
  if (frame) {
    if (is.null(frame.type) || frame.type == "convex") {
      if (is.null(frame.colour) || !(frame.colour %in% 
                                     colnames(plot.data))) {
        hulls <- plot.data[grDevices::chull(plot.data[, 
                                                      1L:2L]), ]
      }
      else {
        hulls <- plot.data %>% dplyr::group_by_(frame.colour) %>% 
          dplyr::do(.[grDevices::chull(.[, 1L:2L]), 
          ])
      }
      mapping <- aes_string(colour = frame.colour, fill = frame.colour)
      p <- p + ggplot2::geom_polygon(data = hulls, mapping = mapping, 
                                     alpha = frame.alpha)
    }
    else if (frame.type %in% c("t", "norm", "euclid")) {
      mapping <- aes_string(colour = frame.colour, fill = frame.colour)
      p <- p + ggplot2::stat_ellipse(mapping = mapping, 
                                     level = frame.level, type = frame.type, geom = "polygon", 
                                     alpha = frame.alpha)
    }
    else {
      stop("frame.type must be convex, t, norm or euclid")
    }
  }
  p <- ggfortify:::post_autoplot(p = p, xlim = xlim, ylim = ylim, log = log, 
                                 main = main, xlab = xlab, ylab = ylab, asp = asp)
  return(p)
}

my_byplot = function (object, data = NULL, scale = 1, x = 1, y = 2, variance_percentage = TRUE, 
                      ...) 
{
  plot.data <- fortify(model = object, data = data)
  plot.data$rownames <- rownames(plot.data)
  if (ggfortify:::is_derived_from(object, "prcomp")) {
    ve <- object$sdev^2/sum(object$sdev^2)
    PC <- paste0("PC", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    loadings.column <- "rotation"
    lam <- object$sdev[c(x, y)]
    lam <- lam * sqrt(nrow(plot.data))
  }
  else if (ggfortify:::is_derived_from(object, "princomp")) {
    ve <- object$sdev^2/sum(object$sdev^2)
    PC <- paste0("Comp.", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    loadings.column <- "loadings"
    lam <- object$sdev[c(x, y)]
    lam <- lam * sqrt(nrow(plot.data))
  }
  else if (ggfortify:::is_derived_from(object, "factanal")) {
    if (is.null(attr(object, "covariance"))) {
      p <- nrow(object$loading)
      ve <- colSums(object$loading^2)/p
    }
    else ve <- NULL
    PC <- paste0("Factor", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    scale <- 0
    loadings.column <- "loadings"
  }
  else if (ggfortify:::is_derived_from(object, "lfda")) {
    ve <- NULL
    PC <- paste0("PC", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    scale <- 0
    loadings.column <- NULL
  }
  else {
    stop(paste0("Unsupported class for autoplot.pca_common: ", 
                class(object)))
  }
  if (scale != 0) {
    lam <- lam^scale
    plot.data[, c(x.column, y.column)] <- t(t(plot.data[, 
                                                        c(x.column, y.column)])/lam)
  }
  plot.columns <- unique(c(x.column, y.column, colnames(plot.data)))
  plot.data <- plot.data[, plot.columns]
  if (!is.null(loadings.column)) {
    loadings.data <- as.data.frame(object[[loadings.column]][, 
    ])
    loadings.data$rownames <- rownames(loadings.data)
    loadings.columns <- unique(c(x.column, y.column, colnames(loadings.data)))
    loadings.data <- loadings.data[, loadings.columns]
  }
  else {
    loadings.data <- NULL
  }
  if (is.null(ve) | !variance_percentage) {
    labs <- PC
  }
  else {
    ve <- ve[c(x, y)]
    labs <- paste0(PC, " (", round(ve * 100, 2), "%)")
  }
  xlab <- labs[1]
  ylab <- labs[2]
  p <- my_ggbiplot(plot.data = plot.data, loadings.data = loadings.data, 
                   xlab = xlab, ylab = ylab, ...)
  return(p)
}
