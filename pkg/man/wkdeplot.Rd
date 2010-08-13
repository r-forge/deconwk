\name{wkde.plot}
\alias{wkde.plot}
\title{
Plot a weighted kernel density estimate
}
\description{
Function to plot a weighted kernel density estimate.
}
\usage{
wkde.plot(y, eval, sigma, h, w, gamma, method = "exact",
          RUG = TRUE, COMPARE = TRUE,
          XLAB = expression(italic(x)), YLAB = "density",
          COL = 1:2, LTY = rep(1, 2), LWD = rep(1, 2), YLIM = NULL)
}
\arguments{
  \item{y}{the observed, contaminated data.}
  \item{eval}{grid on which the deconvolution density estimate be
    calculated.}
  \item{sigma}{the standard deviation of the contaminating (normal)
    distribution.}
  \item{h}{the bandwidth to be used for the weighted kernel density
    estimate}
  \item{w}{the weights to be used.}
  \item{gamma}{the regularisation parameter to be used; either a scalar
    for methods \code{"exact"} and \code{"svm"}, or a vector of values
    from which a suitable value is selected via \eqn{K}-fold
    cross-validation for methods \code{"exact.cv"} and \code{"svm.cv"}}.
  \item{method}{method to be used to solve the quadratic programming
    problem involved in calculating the weights; if \code{"exact"} or
    \code{"exact.cv"} then \code{\link{solveqp}} is used, otherwise
    \code{\link[kernlab]{ipop}} from the \code{kernlab} package is used.} 
  \item{RUG}{logical; if \code{TRUE} a rug is added to the plot
    indicating the location of the observed value and the weight
    attached to each observation.}
  \item{COMPARE}{logical; if \code{TRUE} a kernel density estimate with
    all weights equal to one is added to the plot.}
  \item{XLAB}{label for the x axis.}
  \item{YLAB}{label for the y axis.}
  \item{COL}{colours in which the density estimates should be drawn.}
  \item{LTY}{line type in which the density estimates should be drawn.}
  \item{LWD}{line width with which the density estimates should be drawn.}
  \item{YLIM}{limits for the y axis.}
}
\details{
  If \code{"eval"} is not specified, it defaults to
\code{seq(min(y)-0.1*sd(y), max(y)+0.1*sd(y), length=100)}.

If \code{"w"} is not specified, it defaults to \code{w.hat(y, sigma=sigma, h=h, gamma=gamma, METHOD=method)}.

If \code{"h"} is not specified, it defaults to \code{bw.SJ(y, method="dpi")}.
}
\value{
  Invisible \code{NULL}.  This function is called for its side effect of
  creating a plot.
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\author{
Martin L Hazelton \email{m.hazelton@massey.ac.nz}

Berwin A Turlach \email{berwin@maths.uwa.edu.au}
}
\seealso{
  \code{\link{wkde}}, \code{\link{w.hat}}
}
\examples{
  set.seed(100719)
  sig <- sqrt(29/40)  # Var(Z)/Var(X) = 0.1
  y <- rden(100, DEN=3, sigma=sig)
  wkde.plot(y, sigma=sig, gamma=2.05)
  wkde.plot(y, sigma=sig, gamma=4.4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
