\name{legend2}
\alias{legend2}
\title{A legend which allows longer lines}
\description{
This function can be used to add legends to plots. This is almost identical to the \code{legend} function, accept it has an extra parameter, \code{seg.len} which allows the user to change the lengths of lines shown in legends.
}
\usage{
legend2(x, y = NULL, legend, fill = NULL, col = par("col"), 
    lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"), 
    box.lwd = par("lwd"), box.lty = par("lty"), pt.bg = NA, cex = 1, 
    pt.cex = cex, pt.lwd = lwd, xjust = 0, yjust = 1, x.intersp = 1, 
    y.intersp = 1, adj = c(0, 0.5), text.width = NULL, text.col = par("col"), 
    merge = do.lines && has.pch, trace = FALSE, plot = TRUE, 
    ncol = 1, horiz = FALSE, title = NULL, inset = 0, seg.len = 2)
}
\arguments{
\item{x, y}{the x and y co-ordinates to be used to position the legend.
  They can be specified by keyword or in any way which is accepted by
  \code{\link{xy.coords}}: See Details.}
\item{legend}{a character or \link{expression} vector.
  of length \eqn{\ge 1}{>= 1} to appear in the legend.  Other
  objects will be coerced by \code{\link{as.graphicsAnnot}}.}
\item{fill}{if specified, this argument will cause boxes filled
  with the specified colors (or shaded in the specified colors)
  to appear beside the legend text.}
\item{col}{the color of points or lines appearing in the legend.}
\item{lty, lwd}{the line types and widths for lines appearing in the
  legend.  One of these two \emph{must} be specified for line drawing.}
\item{pch}{the plotting symbols appearing in the legend, either as
  vector of 1-character strings, or one (multi character)
  string.  \emph{Must} be specified for symbol drawing.}
\item{angle}{angle of shading lines.}
\item{density}{the density of shading lines, if numeric and
  positive. If \code{NULL} or negative or \code{NA} color filling
  is assumed.}
\item{bty}{the type of box to be drawn around the legend.  The allowed
  values are \code{"o"} (the default) and \code{"n"}.}
\item{bg}{the background color for the legend box.  (Note that this is
  only used if \code{bty != "n"}.)}
\item{box.lty, box.lwd}{the line type and width for the legend box.}
\item{pt.bg}{the background color for the \code{\link{points}},
  corresponding to its argument \code{bg}.}
\item{cex}{character expansion factor \bold{relative} to current
  \code{par("cex")}.}
\item{pt.cex}{expansion factor(s) for the points.}
\item{pt.lwd}{line width for the points, defaults to the one for
  lines, or if that is not set, to \code{par("lwd")}.}
\item{xjust}{how the legend is to be justified relative to the legend
  x location.  A value of 0 means left justified, 0.5 means centered
  and 1 means right justified.}
\item{yjust}{the same as \code{xjust} for the legend y location.}
\item{x.intersp}{character interspacing factor for horizontal (x) spacing.}
\item{y.intersp}{the same for vertical (y) line distances.}
\item{adj}{numeric of length 1 or 2; the string adjustment for legend
  text.  Useful for y-adjustment when \code{labels} are
  \link{plotmath} expressions.}
\item{text.width}{the width of the legend text in x (\code{"user"})
  coordinates.  (Should be positive even for a reversed x axis.)
  Defaults to the proper value computed by \code{\link{strwidth}(legend)}.}
\item{text.col}{the color used for the legend text.}
\item{merge}{logical; if \code{TRUE}, \dQuote{merge} points and lines but
  not filled boxes.  Defaults to \code{TRUE} if there are points and lines.}
\item{trace}{logical; if \code{TRUE}, shows how \code{legend} does all
  its magical computations.}
\item{plot}{logical.  If \code{FALSE}, nothing is plotted but the
  sizes are returned.}
\item{ncol}{the number of columns in which to set the legend items
  (default is 1, a vertical legend).}
\item{horiz}{logical; if \code{TRUE}, set the legend horizontally
  rather than vertically (specifying \code{horiz} overrides the \code{ncol}
  specification).}
\item{title}{a character string or length-one expression giving a
  title to be placed at the top of the legend.  Other objects will be
  coerced by \code{\link{as.graphicsAnnot}}.}
\item{inset}{inset distance(s) from the margins as a fraction of the
  plot region when legend is placed by keyword.}
\item{seg.len}{ numeric specifying length of lines in legend. }
}
\details{
  Arguments \code{x, y, legend} are interpreted in a non-standard way to
  allow the coordinates to be specified \emph{via} one or two arguments.
  If \code{legend} is missing and \code{y} is not numeric, it is assumed
  that the second argument is intended to be \code{legend} and that the
  first argument specifies the coordinates.

  The coordinates can be specified in any way which is accepted by
  \code{\link{xy.coords}}.  If this gives the coordinates of one point,
  it is used as the top-left coordinate of the rectangle containing the
  legend.  If it gives the coordinates of two points, these specify
  opposite corners of the rectangle (either pair of corners, in any
  order).

  The location may also be specified by setting \code{x} to a single
  keyword from the list \code{"bottomright"}, \code{"bottom"},
  \code{"bottomleft"}, \code{"left"}, \code{"topleft"},
  \code{"top"}, \code{"topright"}, \code{"right"} and
  \code{"center"}. This places the legend on the inside of the plot
  frame at the given location. Partial argument matching is used.  The
  optional \code{inset} argument specifies how far the legend is inset
  from the plot margins.  If a single value is given, it is used for
  both margins; if two values are given, the first is used for \code{x}-
  distance, the second for \code{y}-distance.

  \dQuote{Attribute} arguments such as \code{col}, \code{pch}, \code{lty},
  etc, are recycled if necessary.  \code{merge} is not.

  Points are drawn \emph{after} lines in order that they can cover the
  line with their background color \code{pt.bg}, if applicable.

  See the examples for how to right-justify labels.
}
\value{
  A list with list components
  \item{rect}{a list with components
    \describe{
      \item{\code{w, h}}{positive numbers giving \bold{w}idth and
	\bold{h}eight of the legend's box.}
      \item{\code{left, top}}{x and y coordinates of upper left corner
	of the box.}
    }
  }
  \item{text}{a list with components
    \describe{
      \item{\code{x, y}}{numeric vectors of length \code{length(legend)},
	giving the x and y coordinates of the legend's text(s).}
    }
  }
  returned invisibly.
}
\author{ Richard Pearson (modified from original \pkg{graphics} package function.) }
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth \& Brooks/Cole.

  Murrell, P. (2005) \emph{R Graphics}. Chapman & Hall/CRC Press.
}
\seealso{
  \code{\link[graphics]{legend}}
}
\examples{
x <- seq(-pi, pi, len = 65)
plot(x, sin(x), type = "l", ylim = c(-1.2, 1.8), col = 3, lty = 2)
points(x, cos(x), pch = 3, col = 4)
lines(x, tan(x), type = "b", lty = 1, pch = 4, col = 6)
title("legend(..., lty = c(2, -1, 1), pch = c(-1,3,4), merge = TRUE)",
      cex.main = 1.1)
legend2(-1, 1.9, c("sin", "cos", "tan"), col = c(3,4,6),
       text.col = "green4", lty = c(2, -1, 1), pch = c(-1, 3, 4),
       merge = TRUE, bg = 'gray90', seg.len=6)
}
\keyword{aplot}
