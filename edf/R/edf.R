edf <- function(x, beta, sigma, reps, FUN, ncores = 1, ...){
  stopifnot(class(FUN) == 'function')
  fun_args <- formalArgs(FUN)
  if( !('x' %in% fun_args) ) stop('FUN must have an argument named x')
  if( !('y' %in% fun_args) ) stop('FUN must have an argument named y')
  
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.numeric(x)) stop("x must be a design matrix with only numeric entries")
  sapply(1:length(as.vector(x)), function(i) stopifnot(is.finite(as.vector(x)[i])))
  
  stopifnot(is.numeric(beta))
  if(!is.vector(beta)) stop("beta must be vector")
  sapply(1:length(beta), function(i) stopifnot(is.finite(beta[i])))
  
  if(ncol(x) != length(beta)) stop("the length of beta must match the number of columns in x")
  
  stopifnot(length(sigma) == 1)
  stopifnot(is.numeric(sigma))
  stopifnot(sigma > 0)
  stopifnot(is.finite(sigma))
  
  stopifnot(length(reps) == 1)
  stopifnot(is.numeric(reps))
  stopifnot(reps > 0)
  stopifnot(is.finite(reps))
  stopifnot(reps == round(reps))
  
  n <- nrow(x)
  p <- ncol(x)

  edf_loop <- function(x, beta, sigma, FUN, ...){
	n <- nrow(x)
	e <- rnorm(n, 0, sigma)
	e_test <- rnorm(n, 0, sigma)
	
    y <- x %*% beta + e
    y_test <- x %*% beta + e_test
        
    FUN_fitted <- FUN(x = x, y = y, ...)
    
    if(! is.numeric(FUN_fitted) ) stop("FUN must return numeric values")
    if(nrow(x) != length(as.vector(FUN_fitted))) stop("FUN must return a vector of fitted values of the same length as y")
    
    out_df <- data.frame(cbind(FUN_fitted, y, y_test))
    colnames(out_df) <- c("yhat", "y", "ytest")
    return(out_df)
  }
  
  fit_list <- mclapply(1:reps, function(i) edf_loop(x = x, beta = beta, sigma = sigma, FUN = FUN, ...) , mc.cores = ncores)
  
  fits <-  do.call("cbind", fit_list)

  fits_yhat <- fits[ , seq(1,3*reps, by = 3)]
  fits_y <- fits[ , seq(2,3*reps, by = 3)]
  fits_ytest <- fits[ , seq(3,3*reps, by = 3)]
  
  #Compute EDF
  rowmeans_y <- apply(fits_y, 1, mean)
  rowmeans_yhat <- apply(fits_yhat, 1, mean)
  
  cov_entries <- (fits_y - rowmeans_y%*%t(rep(1,reps)))*(fits_yhat - rowmeans_yhat%*%t(rep(1,reps)))
  
  cov_yyhat <- apply(cov_entries, 1, function(x) sum(x)/(length(x) - 1))
  
  effective_df <- sum(cov_yyhat)/(sigma^2)
   
  y_vec <- unlist(fits_y, use.names = F)
  yhat_vec <- unlist(fits_yhat, use.names = F)
  ytest_vec <- unlist(fits_ytest, use.names = F)
  
  #Compute MSE
  rss <- mean( (y_vec - yhat_vec)^2  )*n
  
  #Compute in-sample test error
  test_rss <- mean( (ytest_vec - yhat_vec)^2 )*n
  
  out_list <- list(edf = effective_df, rss = rss, test_rss = test_rss)
  return(out_list)
}