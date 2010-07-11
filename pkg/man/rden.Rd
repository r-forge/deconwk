\name{rden}
\alias{rden}
\title{
  Simulates from specific (contaminated) distributions
}
\description{
  Simulates from the distributions considered in Hazelton and Turlach
  (2009); details of the distributions (all Gaussian mixtures) are given
  on pages 221--222.
}
\usage{
rden(N, DEN = 1, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{number of observations to be simulated; Should be a single
    number.}
  \item{DEN}{density to simulate from; possible values are 1, 2, 3 and 4
    corresponding to the densities described in the paper.}
  \item{sigma}{the standardard deviation of the contaminating
    measurement error.}
}
\value{
  A vector with the generated random variates.
}
\details{
  The generated random variates are from \eqn{X+Z} where the
  distribution of \eqn{X} is determined by the argument \code{DEN} and
  \eqn{Z} has a normal distribution with mean zero and standard
  deviation \code{sigma}; \eqn{X} and \eqn{Z} are independent.
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\author{
Martin L Hazelton \email{m.hazelton@massey.ac.nz}
}
\keyword{distribution}
