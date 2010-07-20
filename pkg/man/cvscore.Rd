\name{cv.score}
\alias{cv.score}
\title{
  Calculate the CV score for determining regularisation parameter
}
\description{
  Evaluates the cross-validation criterion (11) of Hazelton and Turlach
  (2009). 
}
\usage{
cv.score(y, sigma, h, gamma,
         METHOD=c("exact", "svm"), K=5, verb=FALSE)
}
\arguments{
  \item{y}{the observed values.}
  \item{sigma}{the standard deviation of the contaminating (normal)
    distribution.}
  \item{h}{the smoothing parameter to be used.}
  \item{gamma}{vector of values from which a suitable value is to be
    selected} 
  \item{METHOD}{method to be used to solve the quadratic programming
    problem involved in calculating the weights; if \code{"exact"} then
    \code{\link{solveqp}} is used, otherwise the routine
    \code{\link[kernlab]{ipop}} from the \code{kernlab} package is used.} 
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
Berwin A Turlach \email{berwin@maths.uwa.edu.au}
}
\seealso{
\code{\link{w.hat}}
}
\examples{
  set.seed(100719)
  sig <- sqrt(29/40)  # Var(Z)/Var(X) = 0.1
  y <- rden(100, DEN=3, sigma=sig)
  h <- bw.SJ(y, method="dpi")
  gamma <- exp(seq(from=0, to=6, length=17))

  save.seed <- .Random.seed
  cv1 <- cv.score(y, sigma=sig, h=h, gamma=gamma, METHOD="exact", verb=TRUE)
  plot(log(gamma), cv1, type="b")
  tmp <- getmin(log(gamma), cv1, which="r")
  abline(v=tmp$xmin)

  .Random.seed <- save.seed
  cv2 <- cv.score(y, sigma=sig, h=h, gamma=gamma, METHOD="svm", verb=TRUE)
  plot(log(gamma), cv2, type="b")
  tmp <- getmin(log(gamma), cv2, which="r")
  abline(v=tmp$xmin)

  .Random.seed <- save.seed
  cv1 <- cv.score(y, sigma=sig, h=h, gamma=gamma, METHOD="exact", K=10, verb=TRUE)
  plot(log(gamma), cv1, type="b")
  tmp <- getmin(log(gamma), cv1, which="r")
  abline(v=tmp$xmin)

  .Random.seed <- save.seed
  cv2 <- cv.score(y, sigma=sig, h=h, gamma=gamma, METHOD="svm", K=10, verb=TRUE)
  plot(log(gamma), cv2, type="b")
  tmp <- getmin(log(gamma), cv2, which="r")
  abline(v=tmp$xmin)
}
\keyword{smooth}
