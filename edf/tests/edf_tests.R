library(edf)
library(MASS)
library(tree)

#Tests for edf with prewritten functions
n <- 150
p <- 10
beta <- c(1,0,0,2,0,0,1.5,0,0,4)
rho <- 0.8
sigma_mat <- (1- rho)*diag(p) + matrix(rho,p,p)
x <- mvrnorm(n, rep(0,p), sigma_mat)
sigma <- 1
B <- 100

## Test edf function with linear model
edf_lm <- edf(x, beta, sigma, B, FUN = lm_fitted, formula = y ~ .)
is.numeric(edf_lm$edf)
is.finite(edf_lm$edf)
length(edf_lm) == 3

## Test edf function with pls
edf_pls <- edf(x, beta, sigma, B, FUN = pls_fitted, formula = y ~ ., components = 2)
is.numeric(edf_pls$edf)
is.finite(edf_pls$edf)
length(edf_pls) == 3

## Test edf function with glmnet
edf_ridge <- edf(x, beta, sigma, B, FUN = glmnet_fitted, lambda = 0.1)
is.numeric(edf_ridge$edf)
is.finite(edf_ridge$edf)
length(edf_ridge) == 3

## Test edf function with stepwise regression
edf_step <- edf(x, beta, sigma, B, FUN = step_fitted, formula = y ~.)
is.numeric(edf_step$edf)
is.finite(edf_step$edf)
length(edf_step) == 3

#Test for edf with a user-defined function
tree_fitted <- function(x,y,formula){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:ncol(x), sep = ""))
  tree1 <- tree(formula, data)
  predict(tree1, data[,-1])
}

edf_tree <- edf(x, beta, sigma, B, FUN = tree_fitted, formula = y ~ .)
is.numeric(edf_tree$edf)
is.finite(edf_tree$edf)
length(edf_tree) == 3

#Test passing in bad arguments

## Make sure edf accepts and wrangles a data frame appropriately
x.df <- as.data.frame(x)

edf_lm <- edf(x.df, beta, sigma, B, FUN = lm_fitted, formula = y ~ .)
is.numeric(edf_lm$edf)
is.finite(edf_lm$edf)
length(edf_lm) == 3

## Bad x (contains factors)
fact1 <- rep(c("a","b","c"), each = 10)
num1 <- rnorm(30)
bad_x <- data.frame(fact1, num1)

tryCatch({edf(x = bad_x, beta = c(1,2), sigma = 1, reps = 10, FUN = lm_fitted, formula = y ~.)},
          error=function(e) {
          if (endsWith(conditionMessage(e), "x must be a design matrix with only numeric entries"))
          print("Caught design matrix error")
          })

## Bad beta
bad_beta1 <- c(rep(1,9), NA) #Bad value 1

tryCatch({edf(x = x, beta = bad_beta1, sigma = 1, reps = 10, FUN = lm_fitted, formula = y ~.)},
          error=function(e) {
            if (grepl("beta", conditionMessage(e)))
            print("Caught beta error")
          })

bad_beta2 <- c(rep(1,9), "a") #Bad value 2

tryCatch({edf(x = x, beta = bad_beta2, sigma = 1, reps = 10, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("beta", conditionMessage(e)))
             print("Caught beta error")
         })

bad_beta3 <- rep(1,9) #Length mismatch

tryCatch({edf(x = x, beta = bad_beta3, sigma = 1, reps = 10, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("beta", conditionMessage(e)))
             print("Caught beta error")
         })


## Bad sigma
bad_sigma1 <- NA

tryCatch({edf(x = x, beta = beta, sigma = bad_sigma1, reps = 10, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("sigma", conditionMessage(e)))
             print("Caught sigma error")
         })

bad_sigma2 <- Inf

tryCatch({edf(x = x, beta = beta, sigma = bad_sigma2, reps = 10, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("sigma", conditionMessage(e)))
             print("Caught sigma error")
         })

bad_sigma3 <- -5

tryCatch({edf(x = x, beta = beta, sigma = bad_sigma3, reps = 10, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("sigma", conditionMessage(e)))
             print("Caught sigma error")
         })

bad_sigma4 <- Inf

tryCatch({edf(x = x, beta = beta, sigma = bad_sigma4, reps = 10, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("sigma", conditionMessage(e)))
             print("Caught sigma error")
         })

## Bad reps
bad_reps1 <- 0

tryCatch({edf(x = x, beta = beta, sigma = sigma, reps = bad_reps1, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("reps", conditionMessage(e)))
             print("Caught reps error")
         })

bad_reps2 <- Inf

tryCatch({edf(x = x, beta = beta, sigma = sigma, reps = bad_reps2, FUN = lm_fitted, formula = y ~.)},
         error=function(e) {
           if (grepl("reps", conditionMessage(e)))
             print("Caught reps error")
         })

## Bad FUN
tryCatch({edf(x = x, beta = beta, sigma = sigma, reps = 10, FUN = "blah")},
         error=function(e) {
           if (grepl("FUN", conditionMessage(e)))
             print("Caught FUN error")
         })

tryCatch({edf(x = x, beta = beta, sigma = sigma, reps = 10, FUN = function(z,p,q) z+p+q )},
         error=function(e) {
           if (grepl("FUN", conditionMessage(e)))
             print("Caught FUN error")
         })

