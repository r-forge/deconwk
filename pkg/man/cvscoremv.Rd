\name{cv.score.mv}
\alias{cv.score.mv}
\title{
  Calculate the CV score for determining regularisation parameter for
  the bivariate case
}
\description{
  Evaluates the cross-validation criterion of Hazelton and Turlach
  (2009) for the bivariate case. 
}
\usage{
cv.score.mv(y, Sigma, H, gamma, K=5, verb=FALSE)
}
\arguments{
  \item{y}{matrix with two columns containing the observed values.}
  \item{Sigma}{the variance-covariance matrix of the contaminating (normal)
    distribution.}
  \item{H}{the matrix of smoothing parameter to be used.}
  \item{gamma}{vector of values from which a suitable value is to be
    selected} 
  \item{K}{number of folds to be used if \code{gamma} is chosen by
    cross-validation; defaults to 5.}
  \item{verb}{logical; if \code{TRUE} some progress report will be
    printed during cross-validation.}}
\value{
  A vector containing the cross-validation criterion evaluated at the
  values given in \code{gamma}.
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
\code{\link{w.hat.mv}}
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

  gamma.ridge <- exp(seq(from=-2, to=3, length=17))

  set.seed(100813)
  cv <- cv.score.mv(SBP.bi, Sigma=Sigma.fram, H=H.fram,
                    gamma=gamma.ridge, verb=TRUE)
  plot(log(gamma.ridge), cv, type="b")
  tmp <- getmin(log(gamma.ridge), cv, which="r")
  abline(v=tmp$xmin)

  set.seed(100813)
  cv <- cv.score.mv(SBP.bi, Sigma=Sigma.fram, H=H.fram,
                    gamma=gamma.ridge, verb=TRUE, K=10)
  plot(log(gamma.ridge), cv, type="b")
  tmp <- getmin(log(gamma.ridge), cv, which="r")
  abline(v=tmp$xmin)
}
\keyword{smooth}
