# True values of the model coefficients 
b0 = 0.5 ## intercept
b1 = 0.25 ## log prevalence ratio for X (conditioning on Z)
b2 = 0.1 ## log prevalence ratio for Z (conditioning on X)

# Simulate data 
set.seed(16) ## be a reproducible queen
z = rbinom(n = 100, size = 1, prob = 0.3) ## Z ~ Bern(p = 0.3)
x = runif(100, min = 0, max = 1 + z) ## X|Z ~ Unif(min = 0, max = 1 + Z)
lambda = exp(b0 + b1 * x + b2 * z) ## mean of the Poisson distribution for Y|X,Z
y = rpois(100, lambda = lambda) ## Y|X,Z ~ Pois(lambda), where lambda is a function of X, Z
dat = data.frame(y, x, z) 

# Function to return score vector for Poisson regression 
## Arguments... 
## beta: a vector of parameters (of dimension = dimension(X, Z))
## Y: a number or string identifying the column in data containing outcome Y
## X: a number or string identifying the column in data containing covariate X
## Z: one or more numbers or strings identifying the column(is) in data containing additional covariates Z. 
### If there are no additional covariates, Z = NULL (the default).
## data: a dataframe containing, at least, columns X and Z
U <- function(beta, Y, X, Z = NULL, data) {
  # transform inputs into matrices
  beta = matrix(data = beta, ncol = 1)
  data = data.matrix(frame = data)
  
  # initialize lambda and score
  l <- exp(beta[1] + data[, c(X, Z)] %*% beta[-1]) ## dimension n x 1
  score <- matrix(data = 0, 
                  nrow = nrow(beta)) ## dimension p x 1 (where p = dim(beta))
  
  # save Y - lambda 
  y_min_l = data[, Y] - l ## dimension n x 1 
  
  #fill first element of score
  score[1] <- sum(y_min_l)
  
  # save a transformed version of Y - lambda
  ## same column, replicated once for each covariate in (X, Z)
  y_min_l_wide = matrix(data = y_min_l, 
                        nrow = length(y_min_l),
                        ncol = length(c(X, Z)), 
                        byrow = FALSE) ## dimension n x (p - 1)

  #fill remaining elements of score
  score[-1] <- colSums(data[, c(X, Z)] * y_min_l_wide)
  
  #return score
  score
}

U(beta = c(0.5, 0.25, 0.1), 
  Y = "y", 
  X = "x", 
  Z = "z", 
  data = dat)
# [,1]
# [1,]  2.139859
# [2,] 20.802934
# [3,] 13.943553

# For comparison, here is your old score function with an added "z" argument and a third element of the score
U_old <- function(beta, x, y, z) {
  #initialize lambda, score, sample size
  l <- exp(beta[1] + beta[2]*x + beta[3] * z)
  score <- rep(0, times = 3)
  
  #fill first element of score
  score[1] <- sum(y - l)
  
  #fill second element of score
  score[2] <- sum(x * (y - l))
  
  #fill third element of score
  score[3] <- sum(z * (y - l))
  
  #return score
  score
}

U_old(beta = c(0.5, 0.25, 0.1), 
  y = dat$y, 
  x = dat$x,
  z = dat$z)
# [1]  2.139859 20.802934 13.943553