###### Simulations with Different Covariance Structures #####
library(edf)
library(MASS)
library(Renvlp)

# I typically use PBS scripts to run these simulations and set the number of cores for parallelization to PBS_NUM_PPN (if it is an environment variable on the system).
# Users with different setups may want to set ncores manually
ppn <- as.numeric(Sys.getenv("PBS_NUM_PPN"))
ncores <- ifelse(is.na(ppn), 1, ppn)

# Predictor envelope fitting function 
envlp_fit <- function(x,y,u){
  n <- nrow(x)
  env1 <- xenv(X = x,Y = y, u = u, asy = T)
  sapply(1:n, function(i) pred.xenv(m = env1, Xnew = x[i,])$value )
}

##### Simulation Settings #####
# Number of Monte Carlo iterations
B <- 400

# Number of outer loop iterations
loop_iter <- 100

# Data generation settings
n <- 150
p <- 50
sigma <- 1

rhos <- c(0, 0.2, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8)
dims <- 1:10

# Matrices for output (rows correspond to different rho values, columns correspond to different dim values)
dim_mean_matrix <- cbind(rhos, matrix(0, length(rhos), length(dims)))
dim_var_matrix <- cbind(rhos, matrix(0, length(rhos), length(dims)))
dim_rss_matrix <- cbind(rhos, matrix(0, length(rhos), length(dims)))
dim_test_rss_matrix <- cbind(rhos, matrix(0, length(rhos), length(dims)))

##### Run Simulations #####
for(r in 1:length(rhos)){
  rho <- rhos[r]
  
  # Sigma generation for compound symmetry simulations
  sigma_mat <- (1- rho)*diag(p) + matrix(rho,p,p)
  
  # Sigma generation for block correlation simulations
  #sigma_mat <- matrix(0,p,p)
  #for(i in 1:floor(p/block_size)){
  #	low_ind <- (i-1)*block_size + 1
  #	high_ind <- i*block_size
  #	sigma_mat[(low_ind:high_ind),(low_ind:high_ind)] <- diag(1 - rho, block_size, block_size) + matrix(rho, block_size, block_size)
  #}
  
  # Sigma generation for AR1 simulations
  #sigma_mat <- matrix(0,p,p)
  #for(i in 1:p){
  #	for(j in 1:p) sigma_mat[i,j] <- rho^abs(i-j)
  #}
  
  edf_dims <- matrix(0, nrow = loop_iter, ncol = length(dims))
  rss_dims <- matrix(0, nrow = loop_iter, ncol = length(dims))
  test_rss_dims <- matrix(0, nrow = loop_iter, ncol = length(dims))
  
  for(i in 1:loop_iter){
    x <- mvrnorm(n, rep(0,p), sigma_mat)
    beta_init <- c(rgamma(5,2,scale = 2), rep(0,p-5))
    
    # Shuffle betas for compound symmetry and AR1 simulations, but not block simulations
    beta <- sample(beta_init)
    
    for(j in dims){
      edf_out <- edf(x, beta, sigma, B, FUN = envlp_fit, ncores = ncores, u = j)
      edf_dims[i,j] <- edf_out$edf
      rss_dims[i,j] <- edf_out$rss
      test_rss_dims[i,j] <- edf_out$test_rss
    }
  }
  
  dim_means <- apply(edf_dims, 2, mean)
  dim_mean_matrix[r,-1] <- dim_means
  
  dim_vars <- apply(edf_dims,2, var)
  dim_var_matrix[r,-1] <- dim_vars
  
  rss_means <- apply(rss_dims, 2, mean)
  dim_rss_matrix[r,-1] <- rss_means
  
  test_rss_means <- apply(test_rss_dims, 2, mean)
  dim_test_rss_matrix[r,-1] <- test_rss_means
}

# Because this script takes so long to run, I had it save the output so I could import it later to create plots for the paper.
saveRDS(dim_mean_matrix, "cs_means.Rda")
saveRDS(dim_var_matrix, "cs_vars.Rda")
saveRDS(dim_rss_matrix, "cs_rss.Rda")
saveRDS(dim_test_rss_matrix, "cs_test_rss.Rda")
