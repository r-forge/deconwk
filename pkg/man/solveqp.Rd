\name{solveqp}
\alias{solveqp}
\title{
  Solves a specific quadratic programming problem
}
\description{
  Solves the quadratic programming problem (9) of Hazelton and Turlach
  via a homotopy algorithm approach as described in Appendix B.
}
\usage{
solveqp(Qmat, bvec)
}
\arguments{
  \item{Qmat}{The matrix \bold{Q} in equation (9a) of Hazelton and
    Turlach (2009)}
  \item{bvec}{The vector \bold{\emph{b}} in equation (9a) of Hazelton and
    Turlach (2009)}
}
\value{
  The vector \bold{\emph{w}} that solves the quadratic problem (9).

  Note, the entries in this vector add to one as the code works with a
  different parameterisation of the weight vector.
}
\references{
Hazelton, M.L. and Turlach, B.A. (2009). Nonparametric density
deconvolution by weighted kernel estimators,  Statistics and Computing
19(3): 217--228. \url{http://dx.doi.org/10.1007/s11222-008-9086-7}.  
}
\author{
Berwin A Turlach \email{Berwin.Turlach@gmail.com}
}
\seealso{
  \code{\link[kernlab]{ipop}}
}
\keyword{optimize}
