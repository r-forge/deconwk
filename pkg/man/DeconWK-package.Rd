\name{DeconWK-package}
\alias{DeconWK-package}
\alias{DeconWK}
\docType{package}
\title{
Deconvolution by Weighted Kernels
}
\description{
This package contains code for density deconvolution using weighted
kernel estimators.  Type \sQuote{citation(\dQuote{DeconWK})} for details
of the implemented methods. 
}
\details{
The main functions are:
\tabular{ll}{
\code{\link{w.hat}}:\tab Calculates the weights for density deconvolution using
weighted kernel estimators\cr 
\code{\link{w.hat.mv}}:\tab Calculates the weights for density
deconvolution using bivariate weighted kernel estimators\cr 
\code{\link{wkde}}:\tab Calculates a weighted kernel density estimate\cr
\code{\link{wkde.2d}}:\tab Calculates a bivariate weighted kernel density estimates\cr
\code{\link{decon.f}}:\tab Calculates a classical deconvolution estimate
}
}
\author{
Authors: Martin L Hazelton and Berwin A Turlach

Maintainer: Berwin A Turlach \email{berwin@maths.uwa.edu.au}
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\keyword{ package }
