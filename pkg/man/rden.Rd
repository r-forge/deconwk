\name{rden}
\alias{rden}
\alias{dden}
\title{
  Specific (contaminated) distributions
}
\description{
  Density functions and random generation from the distributions
  considered in Hazelton and Turlach (2009); details of the
  distributions (all Gaussian mixtures) are given on pages 221--222.
}
\usage{
dden(eval, DEN=1, sigma=0)
rden(N, DEN = 1, sigma=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eval}{vector of quantiles.}
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
\examples{
  ##
  ## Figure 1 from paper
  ##
  opar <- par(mfrow=c(2,2))
  eval1 <- seq(-4,4,length=200)
  eval2 <- eval1
  eval3 <- seq(-8,7,length=300)
  eval4 <- seq(-2,30,length=320)

  plot(eval1,dden(eval1,DEN=1),type="l",xlab="",ylab="density")
  lines(eval1,dden(eval1,DEN=1,sigma=0.5),lty=3,lwd=1.25)
  lines(eval1,dden(eval1,DEN=1,sigma=0.5/sqrt(2.5)),lty=2,lwd=1.25)
  lines(eval1,dden(eval1,DEN=1,sigma=0.5*sqrt(2)),lty=4,lwd=1.25)
  title("Density 1")
  plot(eval2,dden(eval2,DEN=2),type="l",xlab="",ylab="density")
  lines(eval2,dden(eval2,DEN=2,sigma=sqrt(51/300)),lty=3,lwd=1.25)
  lines(eval2,dden(eval2,DEN=2,sigma=sqrt(51/300)/sqrt(2.5)),lty=2,lwd=1.25)
  lines(eval2,dden(eval2,DEN=2,sigma=sqrt(51/300)*sqrt(2)),lty=4,lwd=1.25)
  title("Density 2")
  plot(eval3,dden(eval3,DEN=3),type="l",xlab="",ylab="density")
  lines(eval3,dden(eval3,DEN=3,sigma=sqrt(1.8125)),lty=3,lwd=1.25)
  lines(eval3,dden(eval3,DEN=3,sigma=sqrt(1.8125)/sqrt(2.5)),lty=2,lwd=1.25)
  lines(eval3,dden(eval3,DEN=3,sigma=sqrt(1.8125)*sqrt(2)),lty=4,lwd=1.25)
  title("Density 3")
  plot(eval4,dden(eval4,DEN=4),type="l",xlab="",ylab="density")
  lines(eval4,dden(eval4,DEN=4,sigma=sqrt(2.516)),lty=3,lwd=1.25)
  lines(eval4,dden(eval4,DEN=4,sigma=sqrt(2.516)/sqrt(2.5)),lty=2,lwd=1.25)
  lines(eval4,dden(eval4,DEN=4,sigma=sqrt(2.516)*sqrt(2)),lty=4,lwd=1.25)
  title("Density 4")
  par(opar)
}
\keyword{distribution}
