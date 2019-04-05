


####geom_ord_ellipse####
##' add confidence ellipse to ordinary plot produced by ggord
##'
##'
##' @title geom_ord_ellipse
##' @param mapping aes mapping
##' @param ellipse_pro confidence value for the ellipse
##' @param fill color to fill the ellipse, NA by default
##' @param ... additional parameters
##' @return ggplot layer
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 layer
##' @importFrom utils modifyList
##' @export
##' @author Guangchuang Yu
##' @references \url{http://lchblogs.netlify.com/post/2017-12-22-r-addconfellipselda/}
geom_ord_ellipse <- function(mapping = NULL, ellipse_pro = 0.97, fill = NA, ...) {
  default_aes <- aes_(color = ~Groups, group = ~Groups)
  if (is.null(mapping)) {
    mapping <- default_aes
  } else {
    mapping <- modifyList(default_aes, mapping)
  }

  layer(
    geom = "polygon",
    stat = StatOrdEllipse,
    mapping = mapping,
    position = 'identity',
    data = NULL,
    params = list(
      ellipse_pro = ellipse_pro,
      fill = fill,
      ...
    )
  )
}

##' @importFrom ggplot2 ggproto
##' @importFrom ggplot2 Stat
##' @importFrom plyr ddply
##' @importFrom grDevices chull
StatOrdEllipse <- ggproto("StatOrdEllipse", Stat,
                          compute_group = function(self, data, scales, params, ellipse_pro) {
                            names(data)[1:2] <- c('one', 'two')
                            theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
                            circle <- cbind(cos(theta), sin(theta))
                            ell <- ddply(data, .(group), function(x) {
                              if(nrow(x) <= 2) {
                                return(NULL)
                              }
                              sigma <- var(cbind(x$one, x$two))
                              mu <- c(mean(x$one), mean(x$two))
                              ed <- sqrt(qchisq(ellipse_pro, df = 2))
                              data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
                            })
                            names(ell)[2:3] <- c('one', 'two')
                            ell <- ddply(ell, .(group), function(x) x[chull(x$one, x$two), ])
                            names(ell) <- c('Groups', 'x', 'y')
                            return(ell)
                          },
                          required_aes = c("x", "y", "group")
)


## . function was from plyr package
. <- function (..., .env = parent.frame()) {
  structure(as.list(match.call()[-1]), env = .env, class = "quoted")
}


####geom_flat_violin####

## source code from https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R

## somewhat hackish solution to:
## https://twitter.com/EamonCaddigan/status/646759751242620928
## based mostly on copy/pasting from ggplot2 geom_violin source:
## https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r


##' flat violin geom
##'
##'
##' @title geom_flat_violin
##' @param mapping aesthetic mapping
##' @param data the data to display in this layer
##' @param stat The statistical transformation to use on the data for this layer, as a string.
##' @param position position adjustment
##' @param trim If ‘TRUE’ (default), trim the tails of the violins to the range of the data. If ‘FALSE’, don't trim the tails.
##' @param scale if "area" (default), all violins have the same area (before trimming the tails). If "count", areas are scaled proportionally to the number of observations. If "width", all violins have the same maximum width.
##' @param show.legend whether show the legend of this layer
##' @param inherit.aes whether inherit aesthetic mapping from `ggplot`
##' @param ... additional parameters
##' @return ggplot layer
##' @importFrom dplyr group_by
##' @export
##' @author David Robinson
##' @examples
##' library(ggplot2)
##' ggplot(diamonds, aes(cut, carat)) +
##'    geom_flat_violin() +
##'    coord_flip()
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @format NULL
#' @usage NULL
#' @importFrom dplyr mutate
#' @importFrom ggplot2 draw_key_polygon
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)

            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)

          },

          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))

            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))

            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])

            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },

          draw_key = draw_key_polygon,

          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),

          required_aes = c("x", "y")
  )

"%||%" <- getFromNamespace("%||%", "ggplot2")
"%>%"  <- getFromNamespace("%>%", "magrittr")





####geom_cake####
##' @importFrom grid gList
##' @importFrom grid rectGrob
##' @importFrom grid polygonGrob
##' @importFrom grid gpar
candleGrob <- function(x, y, colour.candle = "orange", colour.fire = "red", vp=NULL) {
  width <- 0.02
  height <- 0.2

  xx = c(x+.005,x-.01,x+.01,x+.03,x+.015,x+0.005)
  yy = c(y+.2,y+.23,y+.26,y+.23,y+.2,y+.2)

  gTree(children = gList(
    rectGrob(x+width/2, y+height/2, width = width, height = height, gp = gpar(fill=colour.candle), vp=vp),
    polygonGrob(xx, yy, gp = gpar(fill = colour.fire), vp=vp)
  ))
}

ellipseGrob <- function(x, y, a, b, gp=gpar(), vp=NULL) {
  t <- seq(0, 2*pi, length.out=100)
  xx <- x + a * cos(t)
  yy <- y + b * sin(t)
  polygonGrob(xx, yy, gp = gp, vp=vp)
}

##' @author Guangchuang Yu
##' @importFrom grid segmentsGrob
cakeGrob <- function(x=.5, y=.5, a=.4, b=.14, A=.44, B=.17, height=.3, gp=gpar(), vp=NULL) {
  gp2 <- gp
  if (!is.null(gp$fill)) {
    gp2$col <- gp2$fill
  }
  gTree(children = gList(
    ellipseGrob(x, y-height, A, B, gp=gp, vp=vp),
    ellipseGrob(x, y-height, a, b, gp=gp, vp=vp),
    rectGrob(x, y-height/2, a*2, height, gp=gp2, vp=vp),
    segmentsGrob(x-a, y-height, x-a, y, gp=gp, vp=vp),
    segmentsGrob(x+a, y-height, x+a, y, gp=gp, vp=vp),
    ellipseGrob(x, y, a, b, gp=gp, vp=vp))
  )
}


##' @importFrom grid gTree
cakeCandleGrob <- function(colour.cake = "pink", colour.candle="orange", colour.fire="red", vp=NULL, name=NULL) {
  grobs <- gList(cakeGrob(x=.5, y=.5, a=.4, b=.14, A=.44, B=.17, height=.3, gp=gpar(fill=colour.cake)),
                 candleGrob(.25,.45, colour.candle, colour.fire),
                 candleGrob(.3,.5, colour.candle, colour.fire),
                 candleGrob(.4, .45,colour.candle, colour.fire),
                 candleGrob(.5,.5, colour.candle, colour.fire),
                 candleGrob(.6, .45, colour.candle, colour.fire),
                 candleGrob(.7, .52, colour.candle, colour.fire)
  )
  gTree(children=grobs, name=name, vp=vp)
}


##' ggplot2 layer of birthday cake
##'
##'
##' @title geom_cake
##' @param mapping aes mapping
##' @param data data
##' @param ... additional parameters
##' @return ggplot2 layer
##' @importFrom ggplot2 layer
##' @export
##' @author Guangchuang Yu
geom_cake <- function(mapping = NULL, data = NULL, ...) {
  layer(
    data = data,
    mapping = mapping,
    geom = GeomCake,
    stat = "identity",
    position = "identity",
    params = list(...),
    check.aes = FALSE
  )
}

##' @importFrom grid viewport
##' @importFrom ggplot2 ggproto
##' @importFrom ggplot2 Geom
##' @importFrom ggplot2 draw_key_blank
##' @importFrom ggplot2 aes
GeomCake <- ggproto("GeomCake", Geom,
                    draw_panel = function(data, panel_scales, coord) {
                      data <- coord$transform(data, panel_scales)

                      grobs <- lapply(1:nrow(data), function(i) {
                        vp <- viewport(x=data$x[i], y=data$y[i],
                                       width=data$size[i], height=data$size[i],
                                       angle = data$angle[i],
                                       just = c("center", "center"),
                                       default.units = "native")
                        cakeCandleGrob(data$colour.cake[i], data$colour.candle[i], data$colour.fire[i], vp=vp, name=i)
                      })
                      class(grobs) <- "gList"
                      ggplot2:::ggname("geom_cake",
                                       gTree(children = grobs))
                    },
                    non_missing_aes = c("x", "y", "size", "colour.cake", "colour.candle", "colour.fire"),
                    default_aes = aes(size=.1, colour.cake="#FF3399", colour.candle = "orange", colour.fire="red", angle=0),
                    draw_key = draw_key_blank
)






























