\name{Effective Degrees of Freedom}
\alias{edf}
\title{Compute effective degrees of freedom}
\description{
  The \code{edf} function computes the effective degrees of freedom for models provided by the user using Monte Carlo approximation with linear data generation.
}
\usage{
edf(x, beta, sigma, reps, FUN, ncores, ...)
}
\arguments{
  \item{x}{the design matrix to use for data generation in each Monte Carlo iteration. The user needs to include an intercept column if they want an intercept in the data generation. If the user passes in a matrix of predictors (i.e. a matrix containing factors) instead of a design matrix (with dummy variables in the place of factors), \code{edf} will return an error.}
  \item{beta}{the true coefficient vector to use for data generation}
  \item{sigma}{the standard deviation of the errors in the data generation}
  \item{reps}{the number of Monte Carlo iterations to run}
  \item{FUN}{a model fitting method which returns fitted values. Must take arguments named \code{x} and \code{y}.}
  \item{ncores}{number of cores for parallelization with mclapply}
  \item{...}{additional arguments to pass to FUN}
}
\value{
  \item{edf}{returns the estimated effective degrees of freedom for the model}
  \item{rss}{RSS for fitted values from the model}
  \item{test_rss}{in-sample test RSS for fitted values from the model}
}
\details{
The \code{edf} function uses Monte Carlo to approximate the effective degrees of freedom for a user-specified model fitting method (passed in as the \code{FUN} argument). Examples of valid \code{FUN} arguments and details about how to write them are given by \code{\link{lm_fitted}}, \code{\link{pls_fitted}}, \code{\link{glmnet_fitted}}, and \code{\link{step_fitted}}. The user also specifies the data generation setup \eqn{Y = X\beta + \epsilon} by providing a design matrix \code{x}, a true coefficient vector \code{beta}, and the standard deviation \code{sigma} of the normally distributed error term. 

The effective degrees of freedom for a model fitting method which returns fitted values \eqn{\hat{Y_1}, \ldots, \hat{Y_n}} are given by
  \deqn{df(\hat{Y}) = \frac{1}{\sigma^2} \sum_{i = 1}^{n} Cov(\hat{Y}_i, Y_i) }
(see Hastie et al. 2009). This is a reasonable generalization of the classical concept of degrees of freedom. A brief derivation yields that the effective degrees of freedom for a standard linear model with p-1 predictors, in which \eqn{\hat{Y} = X(X^TX)^{-1}X^TY}, are \eqn{df(\hat{Y}) = tr(X(X^TX)^{-1}X^T) = p}.
Hastie et al. further justify the effective degrees of freedom by showing that they are equal to the expected optimism of the training error rate and so give a reasonable measure of how much a model might be overfit to the training data.

The \code{edf} function approximates the effective degrees of freedom by computing the sample covariance of the replicates of \eqn{Y_i} and \eqn{\hat{Y_i}} for each \eqn{i} over \code{reps} Monte Carlo iterations.

The \code{edf} function uses \link{mclapply} from the \link{parallel} package in order to reduce computation time. On Windows systems, which do not support forking, \link{mclapply} simply calls \link{lapply}.
}
\references{
Hastie, T., Tibshirani, R. & Friedman, J. (2009) \emph{The Elements of Statistical Learning}. 2nd Ed., New York, NY: Springer-Verlag.
}
\examples{
  library(MASS)
  
  n <- 150
  p <- 10
  beta <- c(1,0,0,2,0,0,1.5,0,0,4)
  rho <- 0.8
  sigma_mat <- (1- rho)*diag(p) + matrix(rho,p,p)
  
  x <- mvrnorm(n, rep(0,p), sigma_mat)
  sigma <- 1
  
  #Number of Monte Carlo iterations
  B <- 100
  
  edf_lm <- edf(x, beta, sigma, B,
                FUN = lm_fitted, formula = y ~ .)
  
  edf_lm$edf
}
\keyword{regression}
