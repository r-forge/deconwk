rden <- function(N, DEN=1, sigma=0){
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

dden <- function(eval, DEN=1, sigma=0){
  switch(DEN,
         f <- dnorm(eval,mean=0,sd=sqrt(1+sigma^2)),
         f <- 2*dnorm(eval,sd=sqrt(1+sigma^2))/3+dnorm(eval,sd=sqrt(0.04+sigma^2))/3,
         f <- (dnorm(eval,mean=-2.5,sd=sqrt(1+sigma^2))+dnorm(eval,mean=2.5,sd=sqrt(1+sigma^2)))/2,
         f <- 0.4*dnorm(eval,mean=5,sd=sqrt(5+sigma^2))+0.6*dnorm(eval,mean=13,sd=sqrt(13+sigma^2))
         )
  f
}
