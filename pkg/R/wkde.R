wkde <- function(y, eval=NA, w=NA, h=NA){
  n <- length(y)
  if (is.na(eval[1])) eval <- seq(min(y)-0.1*sd(y),max(y)+0.1*sd(y),length=100)
  if (is.na(w[1])) w <- rep(1,n)
  if (is.na(h)) h <- bw.SJ(y, method="dpi")
  f.hat <- dnorm(outer(y, eval, "-"), sd=h)
  f.hat <- sweep(f.hat, 1, w, "*")
  f.hat <- apply(f.hat, 2, sum)/n
  cbind(x=eval, y=f.hat)
}
