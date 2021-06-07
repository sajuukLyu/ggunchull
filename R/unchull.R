
#' smoothUnchull
#'
#' Get a hull of input set of points. The hull can be non-convex.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param nbin number of points used to shape the hull, default 100
#' @param nsm number of points used to perform convolution, should less than \code{nbin}, default 10
#' @param addsm number of additional times of revolution performed, default 1
#' @param qval quantile of each sector, used to determine the edge of the hull, should less than 1, default 0.95
#' @param sfac expansion size factor, larger value means bigger hull, default 1.5
#'
#' @return a list of coordinates of the hull
#'
#' @export
smoothUnchull <- function(x, y, nbin = 100, nsm = 10, addsm = 1, qval = 0.95, sfac = 1.5) {

  # convert to polar
  polarDim <- cart2polar(x, y)

  # remove illegal points (both dx and dy equal 0)
  illPoint <- polarDim$t == -100
  polarDim$r <- polarDim$r[!illPoint]; polarDim$t <- polarDim$t[!illPoint]

  # find hull around points
  gridNum <- ceiling(nbin * 0.5 * polarDim$t / pi); gridNum[gridNum == 0] <- 1
  upVal <- tapply(polarDim$r, gridNum, function(x) quantile(x, qval))

  # perform convolution
  fullVal <- rep(0, nbin)
  fullVal[as.numeric(names(upVal))] <- upVal
  smVal <- kernapply(fullVal[c((nbin - nsm + 1):nbin, 1:nbin, 1:nsm)], kernel("daniell", nsm))
  if(addsm > 0) {
    for(i in 1:addsm) {
      smVal <- kernapply(smVal[c((nbin - nsm + 1):nbin, 1:nbin, 1:nsm)], kernel("daniell", nsm))
    }
  }

  # convert back
  polar2cart(smVal * sfac, 2 * (1:nbin)/nbin, median(x), median(y))
}

#' cart2polar
#'
#' convert Cartesian coordinates to polar coordinates
#'
#' @param x the x coordinate in Cartesian
#' @param y the y coordinate in Cartesian
#' @param xm the x coordinate of origin point in polar
#' @param ym the y coordinate of origin point in polar
#'
#' @return a list containing r: distance to origin point and t: angle in polar coordinate system
#'
#' @export
cart2polar <- function(x, y, xm = median(x), ym = median(y)) {
  stopifnot(length(x) == length(y))

  dx <- x - xm; dy <- y - ym

  r <- sqrt(dx^2 + dy^2)
  t <- rep(-100, length(x))
  t[dx == 0 & dy > 0] <- pi/2
  t[dx == 0 & dy < 0] <- 3*pi/2
  t[dy == 0 & dx > 0] <- 0
  t[dy == 0 & dx < 0] <- pi

  ind <- which(dx > 0); t[ind] <- atan(dy[ind]/dx[ind])
  ind <- which(dx < 0 & dy > 0); t[ind] <- acos(dx[ind]/r[ind])
  ind <- which(dx < 0 & dy < 0); t[ind] <- pi - asin(dy[ind]/r[ind])
  ind <- which(t < 0); t[ind] <- t[ind] + 2*pi

  list(r = r, t = t)
}

#' polar2cart
#'
#' convert polar coordinates to Cartesian coordinates
#'
#' @param r the distance to origin point in polar
#' @param t the angle in polar
#' @param xm the x coordinate of origin point in polar
#' @param ym the y coordinate of origin point in polar
#'
#' @return a list containing x and y in Cartesian coordinate system
#'
#' @export
polar2cart <- function(r, t, xm, ym) {
  stopifnot(length(r) == length(t))

  x <- r * cospi(t) + xm
  y <- r * sinpi(t) + ym

  list(x = x, y = y)
}

#' stat_unchull
#'
#' A ggplot function used to get a hull of input set of points. The hull can be non-convex.
#'
#' @param nbin number of points used to shape the hull, default 100
#' @param nsm number of points used to perform convolution, should less than \code{nbin}, default 10
#' @param addsm number of additional times of revolution performed, default 1
#' @param qval quantile of each sector, used to determine the edge of the hull, should less than 1, default 0.95
#' @param sfac expansion size factor, larger value means bigger hull, default 1.5
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#'
#' @import ggplot2
#'
#' @export
stat_unchull <- function(
  mapping = NULL, data = NULL, geom = "polygon",
  position = "identity", na.rm = F, show.legend = NA,
  inherit.aes = T, nbin = 100, nsm = 10, addsm = 1, qval = 0.95, sfac = 1.5, ...) {

  layer(
    stat = StatUnchull, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, nbin = nbin, nsm = nsm, addsm = addsm, qval = qval, sfac = sfac, ...)
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @import ggplot2
#' @export
StatUnchull <- ggproto(
  "StatUnchull", Stat,
  required_aes = c("x", "y"),

  compute_group = function(
    data, scales, params, nbin = 100, nsm = 10, addsm = 1, qval = 0.95, sfac = 1.5) {

    chullData <- smoothUnchull(
      data$x, data$y, nbin = nbin, nsm = nsm, addsm = addsm, qval = qval, sfac = sfac)

    data.frame(x = chullData$x, y = chullData$y)
  }
)


