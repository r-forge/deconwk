w.hat <- function(y, sigma, h, gamma, K=5, METHOD="exact"){
  if (missing(h))
    h <- bw.SJ(y,method="dpi")

  n <- length(y)

  ff <- paste(".w.hat.", METHOD, sep="")
  ffargs <- list(y=y, sigma=sigma, h=h, gamma=gamma)
  if(!missing(K)) ffargs <- c(ffargs, K=K)
  w.hat <- do.call(ff, ffargs)
  w.hat
}

.w.hat.exact <- function(y, sigma, h, gamma){
  n <- length(y)
  Y <- outer(y,y,"-")

  lambda <- sqrt(h^2+sigma^2)
  mu <- sqrt(lambda^2+h^2)

  Qmat <- dnorm(Y,sd=sqrt(2)*lambda)
  bvec <- rowMeans(dnorm(Y,sd=mu))

  if(!missing(gamma))
    diag(Qmat) <- diag(Qmat)+gamma/n
  
  w.hat <- solveqp(Qmat, bvec)

  w.hat*n
}

.w.hat.exact.cv <- function(y, sigma, h, gamma, K=5){
  n <- length(y)
  res <- .cv.score1(y, h=h, sigma=sigma, gamrng=gamma, K=K)

  tt <- getmin(log(gamma), res, which="r")
  gamma <- exp(tt$xmin)
  res <- .w.hat.exact(y, sigma=sigma, h=h, gamma=gamma)
  attr(res, "gamma") <- gamma
  res
}

.w.hat.svm1 <- function(y, sigma, h, gamma){
  n <- length(y)
  Y <- outer(y,y,"-")

  lambda <- sqrt(h^2+sigma^2)
  mu <- sqrt(lambda^2+h^2)

  Qmat <- dnorm(Y,sd=sqrt(2)*lambda)
  bvec <- rowMeans(dnorm(Y,sd=mu))

  if(!missing(gamma))
    diag(Qmat) <- diag(Qmat)+gamma/n

  w.hat <- ipop(c=-bvec, H=Qmat,
                A=rep(1,n), b= 1, r=0, l=rep(0,n), u=rep(1,n))
  
  w.hat <- w.hat@primal
  w.hat*n
}

.w.hat.svm1.cv <- function(y, sigma, h, gamma, K=5){
  n <- length(y)
  res <- .cv.score1(y, h=h, sigma=sigma, gamrng=gamma, K=K)

  tt <- getmin(log(gamma), res, which="r")
  gamma <- exp(tt$xmin)
  res <- .w.hat.svm1(y, sigma=sigma, h=h, gamma=gamma)
  attr(res, "gamma") <- gamma
  res
}

.cv.score1 <- function(y, h, sigma, gamrng, K=5){
  n <- length(y)

  ind <- rep(1:K, length=n)
  res <- 0*gamrng

  for(grp in 1:K){
    ii <- which(ind==grp)
    yin <- y[ii]
    yout <- y[-ii]
    Y <- outer(yin, yout, "-")
    M1 <- dnorm(Y, sd=sqrt(h*h+sigma*sigma))
    for(i in seq(along=gamrng)){
      w <- w.hat(yout, h=h, sigma=sigma, gamma=gamrng[i], METHOD="exact")
      res[i] <-  res[i] + mean(log(M1%*%w))
    }
  }

  -res
}
