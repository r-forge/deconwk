\name{getmin}
\alias{getmin}
\title{
  Approximates the minimum of a function given on a grid.
}
\description{
  Approximates the minimum of a function given on a grid.  Quadratic
  approximation around the point where the mimimal function value is
  observed is used (if that point is in the interior).
}
\usage{
getmin(x, y, which="global", count.minima=FALSE, verbose=TRUE)
}
\arguments{
  \item{x}{
    Vector with the x-values at which the function is observed.  Should be
    sorted. 
  }
  \item{y}{
    Vector with the function values.
  }
  \item{which}{
    Defines which miminum we want to find.
    Possible values are "global" for the global minimum, "left" for the
    left-most local minimum and "right" for the right-most local minimum.
    Abbreviations ("g", "r", "gl", etc.) may be used.
  }
  \item{count.minima}{
    If TRUE, the number of local minima in the observed function values is
    returned. 
  }
  \item{verbose}{
    If TRUE, the routine will give a warning if any exceptions occur.
  }
}
\value{
  A list with the following elements is returned:

  \item{xmin}{
    The x-coordinate of the minimum.
  }
  \item{ymin}{
    The approximate value of the function at the minimum.
  }
  \item{nmin}{
    The number of local minima in the y-vector (if requested, otherwise 0).
  }
  \item{excep}{
    Indicates whether an exeption has occurred:
    -1 if the minimum was found at the left end,
    1 if the minimum was found at the right end,
    5 if the minimum was in the middle but the quadratic fit yielded a
    location of the minimum which was outside of the interval defined by
    the three points used for the quadratic fit and
    0 in all other cases.
  }
}
\note{
  The vector x must be sorted.
}
\examples{
        x <- -100:100/50
        y <- x*x
        getmin(x,y)
}
\keyword{optimize}
