#Functions for edf
lm_fitted <- function(x, y, formula){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:ncol(x), sep = ""))
  lm1 <- lm(formula, data)
  lm1$fitted.values
}

glmnet_fitted <- function(x, y, alpha = 0, lambda){
  glmnet1 <- glmnet(x, y, alpha = alpha, lambda = lambda, family = "gaussian")
  predict(glmnet1, newx = x)
}

pls_fitted <- function(x , y, formula, components){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:ncol(x), sep = ""))
  pls1 <- mvr(formula, data = data, center = T, scale = T,
              method = "kernelpls", validation = "none", ncomp = components)
  pls1$fitted.values[,1,components]
}

step_fitted <- function(x, y, formula, k = 2){
  data <- as.data.frame(cbind(y,x))
  colnames(data) <- c("y", paste("x", 1:ncol(x), sep = ""))
  lm1 <- lm(formula, data)
  step1 <- step(lm1, k = k, trace = 0)
  step1$fitted.values
}
