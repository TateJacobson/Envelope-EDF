###### Simulations using Eigenvectors of Cov(X) #####
library(edf)
library(MASS)
library(Renvlp)

# I typically use PBS scripts to run these simulations and set the number of cores for parallelization to PBS_NUM_PPN (if it is an environment variable on the system).
# Users with different setups may want to set ncores manually
ppn <- as.numeric(Sys.getenv("PBS_NUM_PPN"))
ncores <- ifelse(is.na(ppn), 1, ppn)

# Predictor Envelope Fitting function 
envlp_fit <- function(x,y,u){
  n <- nrow(x)
  env1 <- xenv(X = x,Y = y, u = u, asy = T)
  sapply(1:n, function(i) pred.xenv(m = env1, Xnew = x[i,])$value )
}

# Function to generate covariance matrix with given eigenvalues and a given first eigenvector 
generate_cov_mat <- function(eigen_values, first_eigen_vec, p){
  v1 <- first_eigen_vec
  v1 <- v1/sqrt(sum(v1^2)) #normalize first eigenvector
  
  # Generate p-1 random pX1 vectors to orthogonalize with v1
  mat <- matrix( rnorm(p*(p-1)), nrow = p, ncol = p-1)
  mat <- cbind(v1, mat)
  
  # Use Gram-Schmidt to orthogonalize matrix
  for(i in 2:(p-1)){
    proj <- matrix(0, nrow = p, ncol = i-1)
    for(j in 1:i-1){
      proj[,j] <-(sum(mat[,i]*mat[,j])/sum(mat[,j]*mat[,j]))*mat[,j]
    }
    mat[,i] <- mat[,i] - apply(proj,1,sum)
    mat[,i] <- mat[,i]/sqrt(sum(mat[,i]*mat[,i]))
  }
  sigma_mat <- mat%*%diag(eigen_values)%*%t(mat)
  out_list <- list(sigma_mat = sigma_mat, eigen_mat =  mat)
  return(out_list)
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

# The rows of the matrix `eigens` are the eigenvalues for the different simulation settings
indices <- c(100, 25, 10, 5, 1)
eigens <- cbind(indices, matrix(1/2, nrow = length(indices), ncol = p-1))
dims <- 1:10

# Matrices for output (rows correspond to different first eigenvalues, columns correspond to different dim values)
dim_mean_matrix <- cbind(indices, matrix(0, length(indices), length(dims)))
dim_var_matrix <- cbind(indices, matrix(0, length(indices), length(dims)))
dim_rss_matrix <- cbind(indices, matrix(0, length(indices), length(dims)))
dim_test_rss_matrix <- cbind(indices, matrix(0, length(indices), length(dims)))

##### Run Simulations #####
for(e in 1:nrow(eigens)){
  eigen <- eigens[e,]
  eigen_vec1 <- c(rep(1,5), rep(0, p-5))
  
  gen_sigma_mat <- generate_cov_mat(eigen, eigen_vec1, p)
  sigma_mat <- gen_sigma_mat$sigma_mat
  eigen_vec1 <- gen_sigma_mat$eigen_mat[,1]
  eigen_vec2 <- gen_sigma_mat$eigen_mat[,2]
  
  edf_dims <- matrix(0, nrow = loop_iter, ncol = length(dims) )
  rss_dims <- matrix(0, nrow = loop_iter, ncol = length(dims) )
  test_rss_dims <- matrix(0, nrow = loop_iter, ncol = length(dims) )
  
  for(i in 1:loop_iter){
    x <- mvrnorm(n, rep(0,p), sigma_mat)
    
    # Alter the following line to test different generation settings
    beta <- eigen_vec1
    
    for(j in dims){
      edf_out <- edf(x, beta, sigma, B, FUN = envlp_fit, ncores = ncores, u = j)
      
      edf_dims[i,j] <- edf_out$edf
      rss_dims[i,j] <- edf_out$rss
      test_rss_dims[i,j] <- edf_out$test_rss
    }
  }
  
  dim_means <- apply(edf_dims, 2, mean)
  dim_mean_matrix[e,-1] <- dim_means
  
  dim_vars <- apply(edf_dims, 2, var)
  dim_var_matrix[e,-1] <- dim_vars
  
  rss_means <- apply(rss_dims, 2, mean)
  dim_rss_matrix[e,-1] <- rss_means
  
  test_rss_means <- apply(test_rss_dims, 2, mean)
  dim_test_rss_matrix[e,-1] <- test_rss_means
}

# Because this script takes so long to run, I had it save the output so I could import it later to create plots for the paper.
saveRDS(dim_mean_matrix, "eigen_means.Rda")
saveRDS(dim_var_matrix, "eigen_vars.Rda")
saveRDS(dim_rss_matrix, "eigen_rss.Rda")
saveRDS(dim_test_rss_matrix,  "eigen_test_rss.Rda")
