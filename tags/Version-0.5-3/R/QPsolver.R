solveqp <- function(Qmat, bvec){

  EPS <- sqrt(.Machine$double.eps)
##  EPS <- 0.01
  
  n <- length(bvec)
  w <- rep(0, n)
  nu <- - max(bvec)
  ii <- which.max(bvec) 
  jj <- ii              
  kappa <- -bvec - nu


  while(TRUE){
    Qii <- Qmat[ii,ii]
    step <- solve(Qii, rep(1,length(ii)))
    
    u <- rep(0,n)
    u[ii] <- step
    u <- Qmat %*% u
    
    tt1 <- kappa/(1-u)
    tt1[jj] <- Inf
    tt1[tt1<0] <- Inf
    tstep1 <- min(tt1)
    inew <- which.min(tt1)

    tt2 <- -w[ii]/step
    tt2[step>=0] <- Inf
    tstep2 <- min(tt2)
    jnew <- which.min(tt2)
    
    if( is.infinite(tstep1) && is.infinite(tstep2) ){
      wnrm <- sum(w)
      stepnrm <- sum(step)
      tstep <- (1-wnrm)/stepnrm
      w[ii] <- w[ii] + tstep * step
      break
    }
    
    if(tstep1 < tstep2){
      tstep <- tstep1
      case1 <- TRUE
    }else{
      tstep <- tstep2
      case1 <- FALSE
    }

    wold <- w
    woldnrm <- sum(w)
    w[ii] <- w[ii] + tstep*step
    wnrm <- sum(w)

    if( wnrm > 1){
      tstep <- tstep*(1-woldnrm)/(wnrm-woldnrm)
      w <- wold
      w[ii] <- w[ii] + tstep*step
      break
    }

    if(case1){
      jj <- ii <- c(ii, inew)
    }else{
      jj <- ii
      ii <- ii[-jnew]      
    }

    Qw <- Qmat %*% w
    nu <- nu + tstep
    kappa <- kappa - tstep*(1-u)

  }
  Qw <- Qmat %*% w
  nu <- nu + tstep
  kappa <- kappa - tstep*(1-u)
  
  if(max(abs(Qw - bvec - kappa - nu*rep(1,n))) > EPS)
    stop("1 optimality criteria are violated")
  if(any(kappa < -EPS))
    stop("2 optimality criteria are violated")
  if(any(w < -EPS))
    stop("3 optimality criteria are violated")
  if(max(abs(w*kappa)) > EPS)
    stop("4 optimality criteria are violated")
  if(abs((sum(w)-1)*nu) > EPS)
    stop("5 optimality criteria are violated")
  w
}
