\name{w.hat}
\alias{w.hat}
\title{
  Calculate weights for deconvolution
}
\description{
  Routine to calculate the weights for deconvolution via weighted kernel
  density estimates.
}
\usage{
w.hat(y, sigma, h, gamma,
      METHOD=c("exact", "exact.cv", "svm", "svm.cv"), K=5, verb=FALSE)
}
\arguments{
  \item{y}{the observed, contaminated data.}
  \item{sigma}{the standard deviation of the contaminating (normal)
    distribution.}
  \item{h}{the bandwidth to be used for the weighted kernel density
    estimate; if missing the bandwidth returned by \code{bw.SJ(y,
      method="dpi")} will be used.}
  \item{gamma}{the regularisation parameter to be used; either a scalar
    for methods \code{"exact"} and \code{"svm"}, or a vector of values
    from which a suitable value is selected via \eqn{K}-fold
    cross-validation for methods \code{"exact.cv"} and \code{"svm.cv"}.}
  \item{METHOD}{method to be used to solve the quadratic programming
    problem involved in calculating the weights; if \code{"exact"} or
    \code{"exact.cv"} then \code{\link{solveqp}} is used, otherwise
    \code{\link[kernlab]{ipop}} from the \code{kernlab} package is used.} 
  \item{K}{number of folds to be used if \code{gamma} is chosen by
    cross-validation; defaults to 5.}
  \item{verb}{logical; if \code{TRUE} some progress report will be
    printed during cross-validation.}
}
\value{
  A vector containing the weights; if \code{gamma} is chosen by
  cross-validation, the selected value is returned as an attribute.
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
  \code{\link{wkde}}
}
\examples{
  set.seed(100719)
  sig <- sqrt(29/40)  # Var(Z)/Var(X) = 0.1
  y <- rden(100, DEN=3, sigma=sig)
  gamma.ridge <- exp(seq(from=0, to=6, length=17))

  save.seed <- .Random.seed
  w1 <- w.hat(y, sigma=sig, gamma=gamma.ridge, METHOD="exact.cv", verb=TRUE)
  plot(y, w1, type="h")
  attributes(w1)
  .Random.seed <- save.seed
  w2 <- w.hat(y, sigma=sig, gamma=gamma.ridge, METHOD="svm.cv", verb=TRUE)
  plot(y, w2, type="h")
  attributes(w2)

  .Random.seed <- save.seed
  w1 <- w.hat(y, sigma=sig, gamma=gamma.ridge, METHOD="exact.cv", K=10, verb=TRUE)
  plot(y, w1, type="h")
  attributes(w1)
  .Random.seed <- save.seed
  w2 <- w.hat(y, sigma=sig, gamma=gamma.ridge, METHOD="svm.cv", K=10, verb=TRUE)
  plot(y, w2, type="h")
  attributes(w2)
}
\keyword{distribution}
\keyword{smooth}
