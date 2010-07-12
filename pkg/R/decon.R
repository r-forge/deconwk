###########################################
### Classical kernel deconvolution code ###
###########################################

decon.f <- function(y,eval=NA,h=NA,sigma){
  n <- length(y)
  if (is.na(eval[1])) eval <- seq(min(y)-sd(y),max(y)+sd(y),length=100)
  m <- length(eval)
  if (is.na(h)) h <- .decon.pi(y=y,sigma=sigma)
  sig.h <- sigma/h
  YE.h <- outer(y, eval, "-")/h
  K <- matrix(0,nrow=n,ncol=m)
  for (i in 1:n){
    for (j in 1:m){
      fun <- function(u,w){ cos(u*w)*(1-u^2)^3*exp((u*sig.h)^2/2) }               
      K[i,j] <- integrate(fun,lower=-1,upper=1,w=YE.h[i,j])$value
    }
  }
  f.hat <- apply(K, 2, mean)/(2*pi*h)
  cbind(x=eval, y=f.hat)
}

.decon.pi <- function(y,sigma){
  n <- length(y)
  sx <- sqrt(var(y)-sigma^2) 
  theta4 <-  1.851247/sx^9
  abias3 <- function(h){
    term1 <- 6*h^2*theta4
    fn4 <- function(u,w){ u^6*(1-u^2)^6*exp((u*sigma/w)^2) }
    term2int <- integrate(fn4,-1,1,w=h)$value
    term2 <- term2int/(2*pi*n*h^7)
    term1+term2
    }
  h3 <- optimize(f=abias3,lower=sx/sqrt(n),upper=2*sx)$minimum
  char.y <- function(u,y){ 
     apply(cos(outer(u,c(outer(y,y,"-")))),1,mean) 
   }
  fn3 <- function(u,h,y){ u^6*char.y(u/h,y)*(1-u^2)^6*exp((u*sigma/h)^2) }
  theta3 <- integrate(fn3,-1,1,h=h3,y=y)$value/(2*pi*h3^7)
  abias2 <- function(h){
    term1 <- 6*h^2*theta3
    fn2 <- function(u,w){ u^4*(1-u^2)^6*exp((u*sigma/w)^2) }
    term2int <- integrate(fn2,-1,1,w=h)$value
    term2 <- term2int/(2*pi*n*h^5)
    term1+term2
    }
  h2 <- optimize(f=abias2,lower=sx/sqrt(n),upper=2*sx)$minimum
  fn <- function(u,h,y){ u^4*char.y(u/h,y)*(1-u^2)^6*exp((u*sigma/h)^2) }
  theta2 <- integrate(fn,-1,1,h=h2,y=y)$value/(2*pi*h2^5)
  amise <- function(h){
    fn5 <- function(u,w){ (1-u^2)^6*exp((u*sigma/w)^2) }
    term1 <- integrate(fn5,-1,1,w=h)$value/(2*pi*n*h)
    term2 <- h^4*9*theta2
    term1 + term2
    }
  h <- optimize(f=amise,lower=sx/sqrt(n),upper=2*sx)$minimum
  h
}
