\name{wkde.2d}
\alias{wkde.2d}
\title{
Calculate a bivariate weighted kernel density estimate
}
\description{
Calculate a bivariate weighted kernel density estimate.
}
\usage{
wkde.2d(y, Eval, w, H)
}
\arguments{
  \item{y}{the observed values; matrix with two columns.}
  \item{Eval}{two dimensional grid (matrix with two columns) on which
    the deconvolution density estimate is to be calculated.}
  \item{w}{the weights to be used.}
  \item{H}{the matrix of smoothing parameter to be used.}
}
\details{
If \code{"w"} is not specified, it defaults to a vector of ones.

If \code{"H"} is not a matrix, it defaults to \code{Hpi(y)}.
}
\value{
  A vector containing the bivariate deconvolution density estimate. 
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}. 
}
\author{
Martin L Hazelton \email{m.hazelton@massey.ac.nz}
}
\seealso{
 \code{\link{w.hat.mv}}, \code{\link{wkde.contour}},  \code{\link[ks]{Hpi}}
}
\examples{
  library(ks)
  Age <- framingham[,2]
  Age.lim.2 <- 56
  SBP1.A <- framingham[Age>=Age.lim.2,3] # SBP, measure 1, Exam 2 
  SBP2.A <- framingham[Age>=Age.lim.2,4] # SBP, measure 2, Exam 2 
  SBP1.B <- framingham[Age>=Age.lim.2,5] # SBP, measure 1, Exam 3 
  SBP2.B <- framingham[Age>=Age.lim.2,6] # SBP, measure 2, Exam 3 
  sigma.fram.A <- sd(SBP1.A-SBP2.A)
  sigma.fram.B <- sd(SBP1.B-SBP2.B)
  Sigma.fram <- diag(c(sigma.fram.A,sigma.fram.B))^2
  SBP.A <- SBP1.A
  SBP.B <- SBP1.B
  SBP.bi <- cbind(SBP.A,SBP.B)
  H.fram <- Hpi(SBP.bi)

  y1.grid <- seq(min(SBP.bi[, 1]) - 0.5 * sd(SBP.bi[, 1]),
                 max(SBP.bi[, 1]) + 0.5 * sd(SBP.bi[, 1]), length = 25)
  y2.grid <- seq(min(SBP.bi[, 2]) - 0.5 * sd(SBP.bi[, 2]),
                 max(SBP.bi[, 2]) + 0.5 * sd(SBP.bi[, 2]), length = 25)
  Eval <- as.matrix(expand.grid(y1.grid, y2.grid))

  w <-  w.hat.mv(SBP.bi, Sigma.fram, H.fram, gamma = 0.4)
  fhat <- wkde.2d(SBP.bi, Eval = Eval, w = w, H = H.fram)
  str(fhat)
}
\keyword{distribution}
\keyword{smooth}

