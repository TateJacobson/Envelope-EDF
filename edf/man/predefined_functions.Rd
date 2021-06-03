\name{EDF Example Functions}
\alias{lm_fitted}
\alias{pls_fitted}
\alias{glmnet_fitted}
\alias{step_fitted}
\title{
Example functions to pass to \code{edf}.
}
\description{
These functions serve as examples of possible \code{FUN} arguments to pass to the \code{edf} function.
}
\usage{
lm_fitted(x, y, formula)
pls_fitted(x , y, formula, components)
glmnet_fitted(x, y, alpha = 0, lambda)
step_fitted(x, y, formula, k = 2)
}
\arguments{
  \item{x}{covariate matrix}
  \item{y}{response vector}
  \item{formula}{formula for the fitted model}
  \item{components}{the number of components to include in your partial least squares model (see \link[pls]{mvr} for details)}
  \item{alpha}{elastic net mixing parameter (see \code{\link[glmnet]{glmnet}} for details)}
  \item{lambda}{elastic net penalty parameter (see \code{\link[glmnet]{glmnet}} for details)}
  \item{k}{parameter for stepwise penalty; \eqn{ k = 2 } for AIC and \eqn{ k = log(n)} for BIC (see \code{\link[stats]{step}} for details))}
}
\details{
Part of the purpose of these examples is to show users how to write their own functions to pass as \code{FUN} to \code{edf}. As such, I've included the function definitions below for users to use as templates. These functions simply fit models and then return fitted values. Because each model fitting method in R returns fitted values in a different way, it falls on the user to define \code{FUN} in this way.  It is critical that \code{x} and \code{y} are included as named arguments in these functions as \code{edf} specifically uses \code{x} and \code{y} as named arguments in \code{FUN}. 

\preformatted{
lm_fitted <- function(x, y, formula){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:p, sep = ""))
  lm1 <- lm(formula, data)
  lm1$fitted.values
}
}

\preformatted{
pls_fitted <- function(x, y, formula, components){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:p, sep = ""))
  pls1 <- mvr(formula, data = data, center = T, scale = T,
              method = "kernelpls", validation = "none", ncomp = components)
  pls1$fitted.values
}
}

\preformatted{
glmnet_fitted <- (x, y, alpha = 0, lambda){
  glmnet1 <- glmnet(x, y, alpha = alpha, lambda = lambda, family = "gaussian")
  predict(glmnet1, newx = x)
}
}

\preformatted{
step_fitted <- function(x, y, formula, k = 2){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:ncol(x), sep = ""))
  lm1 <- lm(formula, data)
  step1 <- step(lm1, direction = direction, k = k)
  step1$fitted.values
}
}

}
\value{
Each function returns fitted values from a model fit on the data (x,y) provided. These fitted values are then used to compute the effective degrees of freedom in the \code{edf} function.
}

\examples{
library(MASS)

#Generate some data
n <- 150
p <- 50
beta <- c(5,1,0,2.5,0,0,0.5,0, rep(0,37), 10^seq(-1,-5, by = -1))
sigma <- matrix(0,p,p)
for(i in 1:45){
  for(j in 1:45) sigma[i,j] <- 0.5^abs(i-j)
}
for(i in 46:50){
  for(j in 46:50) sigma[i,j] <- 0.1^abs(i-j)
}
set.seed(8053)
x <- mvrnorm(n, rep(0,p), .5*sigma)
e <- rnorm(n, 0, 1)
y <- x \%*\% beta + e

lm_fitted(x,y, formula = y ~ .)
pls_fitted(x,y, formula = y ~ ., components = 2)
glmnet_fitted(x,y,alpha = 0, lambda = 0.01)
step_fitted(x,y, formula = y ~ .)

}
\keyword{models}
\keyword{regression}