rden <- function(N, DEN=1, sigma){
  switch(DEN,
         y <- rnorm(N,sd=sqrt(1+sigma^2)),
         y <- rnorm(N,sd=sqrt( sample(c(0.04,1), size=N,replace=T,prob=c(1,2)/3) + sigma^2)),
         y <- rnorm(N,mean=sample(c(-2.5,2.5),size=N,replace=T),sd=sqrt(1+sigma^2)),
         {
           comp <- sample(1:2,size=N,replace=T,prob=c(.4,.6))
           mu <- c(5,13)
           y <- rnorm(N,mean=mu[comp],sd=sqrt(mu[comp]+sigma^2))
         }
         )
  y
}
