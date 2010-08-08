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
  set.seed(100719)
  sig <- sqrt(29/40)  # Var(Z)/Var(X) = 0.1
  y <- rden(100, DEN=3, sigma=sig)
  h <- bw.SJ(y, method="dpi")
  gamma.ridge <- exp(seq(from=0, to=6, length=17))

  save.seed <- .Random.seed
  cv1 <- cv.score(y, sigma=sig, h=h, gamma=gamma.ridge, METHOD="exact", verb=TRUE)
  plot(log(gamma.ridge), cv1, type="b")
  tmp <- getmin(log(gamma.ridge), cv1, which="r")
  abline(v=tmp$xmin)

  .Random.seed <- save.seed
  cv2 <- cv.score(y, sigma=sig, h=h, gamma=gamma.ridge, METHOD="svm", verb=TRUE)
  plot(log(gamma.ridge), cv2, type="b")
  tmp <- getmin(log(gamma.ridge), cv2, which="r")
  abline(v=tmp$xmin)

  .Random.seed <- save.seed
  cv1 <- cv.score(y, sigma=sig, h=h, gamma=gamma.ridge, METHOD="exact", K=10, verb=TRUE)
  plot(log(gamma.ridge), cv1, type="b")
  tmp <- getmin(log(gamma.ridge), cv1, which="r")
  abline(v=tmp$xmin)

  .Random.seed <- save.seed
  cv2 <- cv.score(y, sigma=sig, h=h, gamma=gamma.ridge, METHOD="svm", K=10, verb=TRUE)
  plot(log(gamma.ridge), cv2, type="b")
  tmp <- getmin(log(gamma.ridge), cv2, which="r")
  abline(v=tmp$xmin)
}
\keyword{smooth}
