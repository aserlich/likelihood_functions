#GET DENOMINATOR
#not sure if this is correct
Dis <- function(lambdas, sigma2) {
  #print(dim(lambdas))
  #print(sigma2)
  #if(class(lambdas) != "matrix") lamdas <- as.matrix(lambdas)
  #obtain theoretical maximum for underdispersion
  #currently uses floor, as suggested in paper
  apply(lambdas,
        MARGIN = 1,
        function(li) sum(seq(0, (-li/(sigma2-1)) + 1))
        )
}

#Get second part of likelihood in two functions
getYi <- function(row, sigma2) {
  #print(row)
  rangeSum <- seq(1, row[1])
  #print(rangeSum)
  quad <- ((sigma2-1) * (rangeSum-1))
  #print(quad)
  out <- sum(log(rep(exp(row[2]), length(quad)) + quad))
  #print(out)
  #print(length(out))
  return(out)
}

secondPart <- function(lambdas, y, sigma2){
  nonzero <- y[y!=0]
  lambdas <- lambdas[y!=0, ,drop=FALSE]
  #print(dim(lambdas))
  #print(length(nonzero))
  allMat <- cbind(nonzero, lambdas)
  out <- apply(allMat, 1, getYi, sigma2)
  #print(out)
}

#Put it all together
gecModel <- function(param, x, y) {
  # preliminaries
  os <- rep(1, nrow(x))
  x <- cbind(os, x)  
  # number of beta estimates equals the number of covariates+constant    
  betas <- param[1:ncol(x)] 
  sigma2 <- param[(ncol(x)+1)]
  lambdas <- x %*% betas
  #print(head(lambdas))
  if (sigma2==1) {
    #niver will happen (there should be a tuning parameter ?)
    Cis <- -exp(lambdas)
  } else {
    if (sigma2 < 1 & sigma2 > 0) {
      denom <- Dis(lambdas, sigma2) 
      Cis <- -exp(lambdas) * log(sigma2) * (sigma2-1)^(-1) - log(denom)
    } else {
      if (sigma2 > 1) {
        #denom <- Dis(lambdas, sigma2) 
        Cis <- -exp(lambdas) * log(sigma2) * (sigma2-1)^(-1) 
      } else
        if (sigma2 <= 0) {
         Cis <- -(abs(lambdas)*100000) #penalty
          warning("There is a penalty")
        } else {
          stop("Error in calculating Cis")
        }
    }
  }
  ll <- -(sum(Cis - y * log(sigma2)) +
             sum(secondPart(lambdas, y, sigma2))
             )
  return(ll)
  }
    #never will empirically happen
    
  # probabilities and penalty function
n <- 10000
x1 <- runif(n)
x2 <- runif(n)
t <- sample(5,n, replace=TRUE)
b0 <- 0
b1 <- 1
b2 <- 2
lambda <- exp(b0 + b1*x1 + b2*x2)
y <- rpois(n, t*lambda)
y2 <- rpois(n, lambda)
y3 <-rnbinom(n, size=2, mu=lambda) # here, size is phi (the dispersion parameter), set beta=3

outPois1 <- optim(c(5,5,5,5), gecModel, hessian = TRUE, x=cbind(x1,x2), y=y2)
outPois2 <- optim(c(5,5,5,.5), gecModel, hessian = TRUE, x=cbind(x1,x2), y=y2)

outOver <- optim(c(5,5,5,2.1), gecModel, hessian = TRUE, x=cbind(x1,x2), y=y3)
