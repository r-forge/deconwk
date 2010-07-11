\name{what}
\alias{w.hat}
\title{
  Calculate weights for deconvolution
}
\description{
  Routine to calculate the weights for deconvolution via weighted kernel
  density estimates.
}
\usage{
w.hat(y, sigma, h, gamma, K = 5, METHOD = "exact")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{METHOD}{
%%     ~~Describe \code{METHOD} here~~
}
}
\value{
  A vector containing the weights.
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\author{
Berwin A Turlach \email{berwin@maths.uwa.edu.au}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{smooth}
