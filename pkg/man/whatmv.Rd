\name{w.hat.mv}
\alias{w.hat.mv}
\title{
  Calculate weights for deconvolution of bivariate KDE
}
\description{
  Routine to calculate the weights for deconvolution via weighted kernel
  density estimates for the bivariate case.
}
\usage{
w.hat.mv(y, Sigma, H, gamma, ...)
}
\arguments{
  \item{y}{a matrix with two columns containing the observed,
    contaminated data.} 
  \item{Sigma}{the variance-covariance matrix of the contaminating
    (normal) distribution.}
  \item{H}{the matrix of smoothing parameters to be used for the
    weighted bivariate kernel density estimate; if missing the bandwidth
    returned by \code{Hpi(y,)} will be used.}
  \item{gamma}{the regularisation parameter to be used; either a scalar
    or a vector of values from which a suitable value is selected via
    \eqn{K}-fold cross-validation.}
  \item{...}{optional parameters passed to the cross-validation routine
    \code{\link{cv.score.mv}} if \eqn{K}-fold cross-validation is used.}
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
  \code{\link{wkde.2d}}, \code{\link{wkde.contour}}, \code{\link[ks]{Hpi}}
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

  w <-  w.hat.mv(SBP.bi, Sigma.fram, H.fram, gamma = 0.4)
  plot(SBP.bi, cex=w)  
}
\keyword{distribution}
\keyword{smooth}
