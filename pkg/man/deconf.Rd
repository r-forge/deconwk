\name{decon.f}
\alias{decon.f}
\title{
  Classical deconvolution density estimate
}
\description{
  Calculates the classical deconvolution density estimate given in
  equation (4) of Hazelton and Turlach (2009).
}
\usage{
decon.f(y, eval = NA, h = NA, sigma)
}
\arguments{
  \item{y}{the observed values.}
  \item{eval}{grid on which the deconvolution density estimate be
    calculated.}
  \item{h}{the smoothing parameter to be used.}
  \item{sigma}{the standard deviation of the contaminating (normal)
    distribution.} 
}
\details{
If \code{"eval"} is not specified, it defaults to
\code{seq(min(y)-sd(y), max(y)+sd(y), length=100)}.

If \code{"h"} is not specified, the plug-in bandwidth selector developed
by Delaigle and Gijbels (2004) is used.
}
\value{
  A matrix with two columns named \code{"x"} and \code{"y"}; the first
  column contains the evaluation grid, \code{"eval"}, and the second
  column the deconvolution density estimate.
}
\references{
Delaigle, A. and Gijbels, I. (2004). Practical bandwidth selection in
deconvolution kernel density estimation.  Computational Statistics &
Data Analysis  45(2): 249--267.
  
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\author{
Martin L Hazelton \email{m.hazelton@massey.ac.nz}
}
\examples{
  set.seed(100712)
  y <- rden(100, DEN=3, sigma=sqrt(29/40)) # Var(Z)/Var(X) = 0.1
  f.hat <- decon.f(y, sigma=sqrt(29/40))
  plot(f.hat, type="l")
}
\keyword{distribution}
\keyword{smooth}
