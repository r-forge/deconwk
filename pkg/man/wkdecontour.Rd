\name{wkde.contour}
\alias{wkde.contour}
\title{
Contour plot of a bivariate weighted kernel density estimate
}
\description{
  Produces a contour plot of a bivariate weighted kernel density estimate
}
\usage{
wkde.contour(y, Sigma, H, w, gamma,
             RUG = TRUE, COMPARE = TRUE, LEVELS = NA,
             XLAB = expression(italic(x)), YLAB = expression(italic(y)),
             DL = FALSE)
}
\arguments{
  \item{y}{the observed values; matrix with two columns.}
  \item{Sigma}{the variance-covariance matrix of the contaminating
    (normal) distribution.}
  \item{H}{the matrix of smoothing parameters to be used for the
    weighted bivariate kernel density estimate; if missing the bandwidth
    returned by \code{Hpi(y,)} will be used.}
  \item{w}{the weights to be used; if missing the weights returned by
    \code{w.hat.mv(y, Sigma, H, gamma = gamma)} will be used.}
  \item{gamma}{the regularisation parameter to be used}
  \item{RUG}{logical; if \code{TRUE} points are added to the plot
    indicating the location of the observed value with the size of the
    points being proportional to the weight attached to each
    observation.} 
  \item{COMPARE}{logical; if \code{TRUE} the contour plot of a kernel
    density estimate with all weights equal to one is added to the plot.}
  \item{LEVELS}{passed to the argument \code{levels} of
    \code{\link[graphics]{contour}}.}
  \item{XLAB}{passed as argument \code{xlab} to
    \code{\link[graphics]{contour}}.}
  \item{YLAB}{passed as argument \code{ylab} to
    \code{\link[graphics]{contour}}.}
  \item{DL}{passed to the argument \code{drawlabs} of
    \code{\link[graphics]{contour}}.}
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

Berwin A Turlach \email{Berwin.Turlach@gmail.com}
}
\seealso{
 \code{\link{wkde.2d}}, \code{\link{w.hat.mv}}, \code{\link[ks]{Hpi}}
}
\examples{
  ##
  ## Figure 7 from paper
  ##
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

  par(mfrow=c(1,2))
  wkde.contour(SBP.bi, Sigma=diag(c(0,0)), H=H.fram,
               RUG=FALSE, COMPARE=FALSE, XLAB="SBP2", YLAB="SBP3",
               LEVELS=seq(5e-5,40e-5,by=10e-5))
  points(SBP.A,SBP.B,pch=19,cex=0.25)

  gamma <- 0.4
  wkde.contour(SBP.bi, Sigma=Sigma.fram, H=H.fram,
               RUG=FALSE, COMPARE=FALSE, XLAB="SBP2", YLAB="SBP3",
               LEVELS=seq(5e-5,40e-5,by=10e-5), gamma=gamma)
  points(SBP.A,SBP.B,pch=19,cex=0.25)
}
\keyword{hplot}
