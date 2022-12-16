#' implementation of uni-variate spline function
SmoothingSpline_Univariate <- function(x, y, lambda){
  N <- NaturalSplineC(x)
  K <- KernelMatrixC(x)
  coef <- solve(t(N)%*%N+lambda*K)%*%t(N)%*%y
  return(coef)
}

#' A implementation of Smoothing Spline function and get the coefficients
#' @param X: Attribute matrix (n*p)
#' @param y: Vector of response
#' @param lambda: Penalty vector on smoothness of each attribute, default to be all 0
#' @param max_iter: Maximum iterations
#' @param tol: Convergence tolerence
#' @return: List of coefficients containing p coefficient vectors
#' @export
SmoothingSpline <- function(X, y, lambda = NA, max_iter = 1000, tol=1e-6){
  Coef <- list()
  
  X <- as.matrix(X)
  if(any(is.na(lambda))){
    lambda <<- rep(0, ncol(X))
  } else{
    if(length(lambda)!=ncol(X)){
      stop("Length of lambda doesn't match")
    }
  }
  # Univariate Case
  if(dim(X)[2]==1){
    Coef[[1]] <- SmoothingSpline_Univariate(X, y, lambda)
    return(Coef)
  } else{
    # Multivariate Case
    n = dim(X)[1]
    p = dim(X)[2]
    # Initialization
    alpha_hat = mean(y)
    f_hat <- matrix(0,nrow=n,ncol=p)
    # Backfitting Algorithm
    flag = 0
    iter = 1
    # First loop to get initialize values of Coef
    for(j in 1:p){
      yj = y - alpha_hat - rowSums(as.matrix(f_hat[,-j]))
      Coef[[j]] = SmoothingSpline_Univariate(X[,j], yj, lambda = lambda[j]) 
      f_hat[,j] <- NaturalSplineC(X[,j])%*%Coef[[j]]
      f_hat[,j] <- f_hat[,j] - mean(f_hat[,j])
    }
    
    while(flag == 0 & iter < max_iter){
      iter = iter + 1
      Coef_old <- Coef
      Coef_diff <- 0
      ## fun
      Coef_diff = calCoefDiffC(p, X, y, alpha_hat, f_hat, lambda, Coef, Coef_old, Coef_diff)
      # Check Convergence
      if(Coef_diff < tol){
        flag <- 1
      }
    }
    if(flag == 0){
      warning("Backfitting algorithm doesn't converge")
    }
    print(paste0("Number of iterations: ",iter))
    return(Coef)
  }
}