\name{framingham}
\alias{framingham}
\docType{data}
\title{
Framingham heart study data
}
\description{
  This is the Framingham data from Carroll et al. (2006)
}
\usage{framingham}
\format{
  A data frame with 1615 observations on the following 14 variables.
  \describe{
    \item{\code{id}}{observation number (1--1615)}
    \item{\code{age}}{Age at exam 2}
    \item{\code{sbp21}}{First systolic blood pressure at exam 2}
    \item{\code{sbp22}}{Second systolic blood pressure at exam 2}
    \item{\code{sbp31}}{First systolic blood pressure at exam 3}
    \item{\code{sbp32}}{Second systolic blood pressure at exam 3}
    \item{\code{smoke}}{Smoking status at exam 1 (1=smoker)}
    \item{\code{cholest2}}{serum cholesterol at exam 2}
    \item{\code{cholest3}}{serum cholesterol at exam 3}
    \item{\code{firstchd}}{First evidence of coronary heart disease
      (CHD) occurring at exam 3 through 6, i.e., within an eight-year
      follow-up period to exam 2 (1=yes)}
  }
}
\details{
  \enumerate{
    \item The data are for \emph{males} only.
    \item The data contain complete records only.
  }
}
\source{
  \url{http://www.stat.tamu.edu/~carroll/eiv.SecondEdition/}

  \url{http://www.stat.tamu.edu/~carroll/eiv.SecondEdition/data.php}
}
\references{
  Carroll, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu,
  C.M. (2006). \emph{Measurement Error in Nonlinear Models: A Modern
    Perspective (2nd ed)}, Chapman & Hall/CRC.
}
\examples{
str(framingham)
}
\keyword{datasets}
