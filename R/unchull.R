#' stat_unchull
#'
#' A ggplot function used to get a hull of input set of points. The hull can be non-convex.
#'
#' @param n parameter for kNN algorithm to defining the point density, default 10
#' @param delta the distance to extend the curve
#' @param th the distance to simplify the curve
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#'
#' @import ggplot2
#' 
#' @name stat_unchull
#' @rdname stat_unchull
#'
#' @export
stat_unchull <- function(
    mapping = NULL, data = NULL, geom = "polygon",
    position = "identity", na.rm = F, show.legend = NA,
    inherit.aes = T, n = 10, delta, th = delta, ...) {
  
  layer(
    stat = StatUnchull, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(n = n, delta = delta, th = th, ...)
  )
}

#' @format NULL
#' @usage NULL
#' @import ggplot2
#' @export
StatUnchull <- ggproto(
  "StatUnchull", Stat,
  required_aes = c("x", "y"),
  
  compute_group = function(data, scales, params, n = 10, delta, th) {
    x <- data$x
    y <- data$y
    
    tmp1 <- step1_denoise(x, y, n)
    tmp2 <- step2_gethull(tmp1$x, tmp1$y, tmp1$nnDist * 2)
    tmp3 <- step3_extend(tmp2$x, tmp2$y, delta)
    tmp4 <- step4_simplify(tmp3$x, tmp3$y, th)
    tmp5 <- step5_smooth(tmp4$x, tmp4$y)
    
    data.frame(x = tmp5$x, y = tmp5$y)
  }
)

#' step1_denoise
#' 
#' Denoising by keeping only the largest cell cluster using dbscan algorithm.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param n parameter for kNN algorithm to defining the point density
#'
#' @import dbscan
#' @import magrittr
#' 
#' @export
step1_denoise <- function(x, y, n) {
  xy <- cbind(x, y)
  
  knnRes <- dbscan::kNN(xy, k = n, sort = T)
  nnDist <- mean(knnRes$dist[, n])
  
  dbscanRes <- dbscan::dbscan(xy, eps = nnDist*2, minPts = n)
  
  maxClus <- dbscanRes$cluster %>% table %>% {names(.)[which.max(.)]}
  xy <- xy[as.character(dbscanRes$cluster) == maxClus, ]
  
  list(x = xy[, 1], y = xy[, 2], nnDist = nnDist)
}

#' step2_gethull
#' 
#' Get the boundary of given cell cluster using Alpha Shape algorithm.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param alpha value of alpha
#'
#' @import alphahull
#'
#' @export
step2_gethull <- function(x, y, alpha) {
  xy <- cbind(x, y)
  
  ashapeRes <- alphahull::ashape(xy, alpha = alpha)
  
  ahull <- ashapeRes$edges[, c("ind1", "ind2")]
  n <- nrow(ahull)
  for(i in 2:n) {
    is_next <- apply(ahull[i:n, , drop = F], 1, function(x) any(ahull[i-1, 2] %in% x))
    next_idx <- which(is_next) + i-1
    
    ahull <- ahull[c(1:(i-1), next_idx, setdiff(i:n, next_idx)), ]
    if(ahull[i, 1] != ahull[i-1, 2]) {
      ahull[i, ] <- ahull[i, 2:1]
    }
  }
  xy <- xy[ahull[, 1], ]
  
  list(x = xy[, 1], y = xy[, 2])
}

#' step3_extend
#' 
#' Extend the boundary of given cell cluster using polygon-clipping.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param delta distance to shift the boundary
#' 
#' @import polyclip
#' 
#' @export
step3_extend <- function(x, y, delta) {
  polyclip::polyoffset(list(x = x, y = y), delta = delta, jointype = "square")[[1]]
}

#' step4_simplify
#' 
#' Simplify the boundary of given cell cluster using Polyline Simplification algorithm.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param th distance to distribute vertexes
#'
#' @export
step4_simplify <- function(x, y, th) {
  xy <- cbind(x, y)
  
  i = 2
  while(i < nrow(xy)) {
    is_tooclose <- (xy[i:nrow(xy), 1] - xy[i-1, 1])^2 + (xy[i:nrow(xy), 2] - xy[i-1, 2])^2 < th^2
    xy <- xy[c(1:(i-1), (i:nrow(xy))[!is_tooclose]), ]
    
    i = i+1
  }
  
  list(x = xy[, 1], y = xy[, 2])
}

#' step5_smooth
#' 
#' Get an X-spline using boundary polygon as control points.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' 
#' @import graphics
#'
#' @export
step5_smooth <- function(x, y) {
  
  tf <- tempfile(fileext=".png")
  png(tf)
  plot.new()
  tmp <- graphics::xspline(x = x, y = y, shape = 1, open = F, draw = F)
  invisible(dev.off())
  unlink(tf)
  
  list(x = tmp$x, y = tmp$y)
}


#' smoothUnchull
#'
#' Get a hull of input set of points. The hull can be non-convex.
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param nbin number of points used to shape the hull, default 100
#' @param nsm number of points used to perform convolution, should less than \code{nbin}, default 10
#' @param addsm number of additional times of convolution performed, default 1
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

#' stat_unchull0
#'
#' A ggplot function used to get a hull of input set of points. The hull can be non-convex.
#'
#' @param nbin number of points used to shape the hull, default 100
#' @param nsm number of points used to perform convolution, should less than \code{nbin}, default 10
#' @param addsm number of additional times of convolution performed, default 1
#' @param qval quantile of each sector, used to determine the edge of the hull, should less than 1, default 0.95
#' @param sfac expansion size factor, larger value means bigger hull, default 1.5
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#'
#' @import ggplot2
#'
#' @export
stat_unchull0 <- function(
  mapping = NULL, data = NULL, geom = "polygon",
  position = "identity", na.rm = F, show.legend = NA,
  inherit.aes = T, nbin = 100, nsm = 10, addsm = 1, qval = 0.95, sfac = 1.5, ...) {

  layer(
    stat = StatUnchull0, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, nbin = nbin, nsm = nsm, addsm = addsm, qval = qval, sfac = sfac, ...)
  )
}

#' @format NULL
#' @usage NULL
#' @import ggplot2
#' @export
StatUnchull0 <- ggproto(
  "StatUnchull", Stat,
  required_aes = c("x", "y"),

  compute_group = function(
    data, scales, params, nbin = 100, nsm = 10, addsm = 1, qval = 0.95, sfac = 1.5) {

    chullData <- smoothUnchull(
      data$x, data$y, nbin = nbin, nsm = nsm, addsm = addsm, qval = qval, sfac = sfac)

    data.frame(x = chullData$x, y = chullData$y)
  }
)


