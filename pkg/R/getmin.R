getmin<-function(x, y, which="global", count.minima=FALSE, verbose=TRUE){

  n    <- length(x)
  if( length(y) != n)
    stop("x and y do not match.")

  flag <- charmatch(which,c("left","global","right"), nomatch=-1)
  if( flag == -1 )
    stop("Parameter `which' must be either `left', `global' or `right'.")
  flag <- flag-2

  xmin <- 0
  ymin <- 0
  if(count.minima)
    nmin <- 1
  else
    nmin <- 0

  res <- .Fortran(.DCWK_getmin,
                  as.double(x), as.double(y), as.integer(n),
                  xmin=as.double(xmin), ymin=as.double(ymin),
                  flag=as.integer(flag), nmin=as.integer(nmin))
  excep <- res$flag

  if(verbose){
    if(excep == -1)
      warning("Minima was found at left end.")
    if(excep == 1)
      warning("Minima was found at right end.")
    if(excep == 5)
      warning("Singularity in quadratic fit.")
  }

  list(xmin=res$xmin, ymin=res$ymin, nmin=res$nmin, excep=excep)
}
