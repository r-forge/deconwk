
w.hat.mv <- function(y, Sigma, H=NA, gamma, ...){
  if (!is.matrix(y)) stop("Use w.hat for univariate estimation")
  if (!is.matrix(H)) H <- Hpi(y)
  n <- dim(y)[1]

  CV <- FALSE
  if(length(gamma) > 1){
    CV <- TRUE
    res <- cv.score.mv(y, Sigma=Sigma, H=H, gamma=gamma, ...)

    tt <- getmin(log(gamma), res, which="r")
    gamma <- exp(tt$xmin)
  }
  
  yy <- y[rep(1:n,n),]-y[rep(1:n,rep(n,n)),]
  Lambda <- H + Sigma
  Omega <- Lambda + H

  Qmat <- matrix(dmvnorm(yy,sigma=2*Lambda),ncol=n)
  if(!missing(gamma))
    diag(Qmat) <- diag(Qmat) + gamma/n
  bvec <- rowMeans(matrix(dmvnorm(yy,sigma=Omega),ncol=n))

  w.hat <- ipop(c=-bvec, H=Qmat,
                A=rep(1,n), b= 1, r=0, l=rep(0,n), u=rep(1,n))
  w.hat <- w.hat@primal*n

  if(CV)
    attr(w.hat, "gamma") <- gamma

  w.hat
}

cv.score.mv <- function(y, Sigma, H, gamma, K=5, verb=FALSE){
  if (!is.matrix(y)) stop("y not multivariate")
  n <- dim(y)[1]

  ind <- sample(rep(1:K, length=n))
  res <- 0*gamma

  for(grp in 1:K){
    if(verb)
      cat("Fold", grp, ":\t")
    ii <- which(ind==grp)
    yin <- y[ii,]
    nin <- length(ii)
    yout <- y[-ii,]
    
    yy <- yin[rep(1:nin, each=n-nin),] - yout[rep(1:(n-nin), nin),]
    M1 <- dmvnorm(yy, sigma=H+Sigma)
    M1 <- matrix(M1, nrow=nin, byrow=TRUE)
    for(i in seq(along=gamma)){
      if(verb)
        cat(i, " ")
      w <- w.hat.mv(yout, H=H, Sigma=Sigma, gamma=gamma[i])
      res[i] <-  res[i] + mean(log(M1%*%w))
    }
  }
  
  -res
}
