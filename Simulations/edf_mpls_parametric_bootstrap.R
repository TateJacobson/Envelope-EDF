###### Parametric Bootstrap Simulation with Minneapolis Schools Data #####
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

##### Simulation Settings #####
# Number of Monte Carlo iterations
B <- 400

# Load Minneapolis Schools data
mpls <- read.table("mpls.txt", header = T, sep = " ")
x <- data.matrix( mpls[, -c(1:4,14)] ) 
y <- as.vector(mpls[, 4])

# Centering and scaling the predictor matrix
col_means <- apply(x, 2, mean)
col_vars <- apply(x, 2, var)
for(i in 1:ncol(x)){
  x[,i] <- (x[,i]-col_means[i])/(sqrt(col_vars[i]))
}

# Determine data generating model with envelope (use column 4 as response)
env_check <- xenv(X = x, Y = y, u = 1, asy = T)
sigma <- sqrt(env_check$SigmaYcX)
beta <- as.vector(env_check$beta)

dims <- 1:5
dim_means <- rep(0, length(dims))

##### Run Simulations #####
for(j in 1:length(dims)){
    dim_means[j] <- edf(x, beta, sigma, B, FUN = envlp_fit, ncores = ncores, u = dims[j])$edf
}

# Because this script takes so long to run, I had it save the output so I could import it later to create plots for the paper.
saveRDS(dim_means, "mpls_bootstrap_means.Rda")
