wkde <- function(y, eval, w, h){
  n <- length(y)
  if (missing(eval))
    eval <- seq(min(y)-0.1*sd(y),max(y)+0.1*sd(y),length=100)
  if (missing(w))
    w <- rep(1,n)
  if (missing(h))
    h <- bw.SJ(y, method="dpi")
  f.hat <- dnorm(outer(y, eval, "-"), sd=h)
  f.hat <- sweep(f.hat, 1, w, "*")
  f.hat <- apply(f.hat, 2, sum)/n
  cbind(x=eval, y=f.hat)
}

wkde.plot <- function(y, eval, sigma, h, w, gamma,
                      method="exact",
                      RUG=TRUE, COMPARE=TRUE,
                      XLAB=expression(italic(x)),
                      YLAB="density",
                      COL=1:2, LTY=rep(1,2), LWD=rep(1,2), YLIM=NULL){
  n <- length(y)
  if(missing(h))
    h <- bw.SJ(y, method="dpi")
  if(missing(eval))
    eval <- seq(min(y)-0.1*sd(y),max(y)+0.1*sd(y),length=100)
  if (missing(w))
    w <- w.hat(y, sigma=sigma, h=h, gamma=gamma, METHOD=method)
  fhat <- wkde(y, eval, w, h)

  plot(fhat,
       xlab=XLAB, ylab=YLAB,
       type="l", col=COL[1], lty=LTY[1], lwd=LWD[1], ylim=YLIM)

  W <- w*0.05/max(w)
  if (RUG) for (i in 1:n) rug(y[i],ticksize=W[i])
  if (COMPARE)
    lines(wkde(y, eval=eval, h=h),
          col=COL[2], lty=LTY[2], lwd=LWD[2])
  invisible(NULL)
}

wkde.2d <- function(y, Eval, w, H){
  n <- dim(y)[1]
  m <- dim(Eval)[1]
  if (missing(w))
    w <- rep(1,n)
  if (!is.matrix(H))
    H <- Hpi(y)
  tmp.index <- rep(1:n,rep(m,n))
  Y <- y[tmp.index,]
  W <- w[tmp.index]
  EVAL <- Eval[rep(1:m,n),]
  YE <- Y-EVAL
  f.hat <- W*dmvnorm(YE,mean=c(0,0),sigma=H)
  f.hat <- c(tapply(f.hat,factor(rep(1:m,n)),mean))
  f.hat
}

wkde.contour <- function(y, Sigma, H, w, gamma,
                         RUG=TRUE, COMPARE=TRUE, LEVELS=NA,
                         XLAB=expression(italic(x)),
                         YLAB=expression(italic(y)), DL=FALSE){
  n <- dim(y)[1]
  if (!is.matrix(H)) H <- Hpi(y)
  y1.grid <- seq(min(y[,1])-0.5*sd(y[,1]),max(y[,1])+0.5*sd(y[,1]),length=25)
  y2.grid <- seq(min(y[,2])-0.5*sd(y[,2]),max(y[,2])+0.5*sd(y[,2]),length=25)
  Eval <- as.matrix(expand.grid(y1.grid,y2.grid))		
  if (missing(w)) w <- w.hat.mv(y,Sigma,H,gamma=gamma)
  fhat <- wkde.2d(y, Eval=Eval, w=w, H=H)
  if (is.na(LEVELS[1])) LEVELS <- pretty(range(fhat),8-4*COMPARE)

  contour(y1.grid, y2.grid, matrix(fhat,ncol=25),
          xlab=XLAB, ylab=YLAB, levels=LEVELS, drawlabels=DL)
  if (RUG) points(y,pch=13,cex=w/2)
  if (COMPARE) {
    contour(y1.grid, y2.grid, matrix(wkde.2d(y, Eval=Eval, H=H),ncol=25),
            col="red", add=TRUE, levels=LEVELS, drawlabels=DL)
    points(y,pch=1,col="red",cex=0.5)
  }
  invisible(NULL)
} 
