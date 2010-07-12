\name{wkde}
\alias{wkde}
\title{
  Weighted kernel density estimate
}
\description{
  Calculates a weighted kernel density estimate as defined by equation
  (5) of Hazelton and Turlach (2009). 
}
\usage{
wkde(y, eval = NA, w = NA, h = NA)
}
\arguments{
  \item{y}{the observed values.}
  \item{eval}{grid on which the deconvolution density estimate be
    calculated.}
  \item{w}{the weights to be used.}
  \item{h}{the smoothing parameter to be used.}
}
\details{
If \code{"eval"} is not specified, it defaults to
\code{seq(min(y)-0.1*sd(y), max(y)+0.1*sd(y), length=100)}.

If \code{"w"} is not specified, it defaults to a vector of ones.

If \code{"h"} is not specified, it defaults to \code{bw.SJ(y, method="dpi")}.
}
\value{
  A matrix with two columns named \code{"x"} and \code{"y"}; the first
  column contains the evaluation grid, \code{"eval"}, and the second
  column the deconvolution density estimate.
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\author{
Martin L Hazelton \email{m.hazelton@massey.ac.nz}
}
\examples{
  set.seed(100712)
  sig <- sqrt(29/40)  # Var(Z)/Var(X) = 0.1
  y <- rden(100, DEN=3, sigma=sig)
  f.hat <- wkde(y)
  plot(f.hat, type="l", ylim=c(0, 0.2))
  w <- w.hat(y, sigma=sig, gamma=1.6) 
  fd.hat <- wkde(y, w=w)
  lines(fd.hat, col="red")
  w <- w.hat(y, sigma=sig, gamma=2.7)
  fd.hat <- wkde(y, w=w)
  lines(fd.hat, col="blue")
}
\keyword{distribution}
\keyword{smooth}
