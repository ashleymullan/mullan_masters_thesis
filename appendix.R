# This function returns the maximum likelihood estimates (MLEs) for 
# the Poisson regression model with covariate misclassification from 
# Mullan et al. (2024+). It calls a helper function, loglik_mat, to 
# evaluate the log likelihood.
#
# Arguments:
# error_formula: a misclassification model formula (or coercible to 
# formula), a formula expression as for other regression models. 
# The response should be the error-free version of the error-prone of 
# the covariate. 
# analysis_formula: a analysis model formula (or coercible to formula), 
# a formula expression as for other regression models. The response 
# should be the Poisson model outcome, and, if needed, the offset can 
# be provided as the offset argument.
# offset: optional, variable name string for the analysis model offset
# that defaults to NULL
# data: a dataset containing at least the variables included above
#
# A list containing a dataframe with final coefficient and standard error 
#estimates for the analysis model and a convergence code are returned

mlePossum = function(error_formula, analysis_formula, offset = NULL, data) {
  ## Extract variable names from user-specified formulas 
  get_Y_name = as.character(as.formula(analysis_formula))[2]
  get_X_name = as.character(as.formula(error_formula))[2]
  
  analysis_covar = unlist(strsplit(x = gsub(pattern = " ", 
                                            replacement = "", 
                                            x = as.character(as.formula(
                                              analysis_formula))[3]), 
                                   split = "+", 
                                   fixed = TRUE))
  error_covar = unlist(strsplit(x = gsub(pattern = " ", 
                                         replacement = "", 
                                         x = as.character(as.formula(
                                           error_formula))[3]), 
                                split = "+", 
                                fixed = TRUE))
  get_Xstar_name = setdiff(error_covar, analysis_covar) 
  
  get_Z_name = intersect(error_covar, analysis_covar) 
  
  ## Fit complete-case models to get initial values 
  ### P(Y|X,Z) 
  if (is.null(offset)) {
    cc_fit = glm(formula = as.formula(analysis_formula), 
                 data = data, 
                 family = poisson)$coefficients
  } else {
    cc_fit = glm(formula = as.formula(paste0(
      paste(get_Y_name, paste(get_X_name, collapse = "+"), sep = "~"), 
      "+offset(log(", offset, "))")), 
                 data = data, 
                 family = poisson)$coefficients
  }
  ### P(X|X*,Z)
  cc_fit = c(cc_fit, 
             glm(formula = as.formula(error_formula), 
                 data = data, 
                 family = binomial)$coefficients)
  
  ## Add queried/non-missing data indicator
  data[, "Q"] = as.numeric(!is.na(data[, get_X_name]))
  
  #Obtain MLEs
  if (length(get_Z_name) > 0) { #case with an extra covariate
    optim_res = optim(fn = loglik_mat, 
                      par = cc_fit, 
                      hessian = TRUE, 
                      method = "BFGS",
                      Y_name = get_Y_name,
                      X_name = get_X_name,
                      Z_name = get_Z_name,
                      Xstar_name = get_Xstar_name,
                      Q_name = "Q",
                      offset_name = offset,
                      data = data)
  } else { #case without an extra covaraite
    optim_res = optim(fn = loglik_mat, 
                      par = cc_fit, 
                      hessian = TRUE, 
                      method = "BFGS",
                      Y_name = get_Y_name,
                      X_name = get_X_name,
                      Xstar_name = get_Xstar_name,
                      offset_name = offset,
                      Q_name = "Q",
                      data = data)
  }
  
  ## Prepare model output to be returned 
  ### Invert the hessian to get estimated standard errors
  num_analysis_covar = length(analysis_covar) + 1
  cov_theta = 
    tryCatch(expr = 
               solve(optim_res$hessian)[1:num_analysis_covar, 
                                        1:num_analysis_covar],
                       error = function(err) { 
                         #return variance/covariance of NAs if Hessian
                         #is not invertible
                         matrix(data = NA, 
                                nrow = num_analysis_covar, 
                                ncol = num_analysis_covar)
                       })
  
  if (optim_res$convergence == 0) { #case where model converges
    res = data.frame(Est = optim_res$par[1:num_analysis_covar],
                     SE = sqrt(diag(cov_theta)))
    rownames(res) = c("(Intercept)", analysis_covar)
  } else { #case where model does not converge
    res = data.frame(Est = rep(NA, num_analysis_covar),
                     SE = rep(NA, num_analysis_covar))
    rownames(res) = c("(Intercept)", analysis_covar)
  }
  
  ## Return model output and convergence
  return(list(coefficients = res, 
              convergence = optim_res$convergence))
}


# This is the helper function. The parameters are extracted from 
# information obtained by the mlePossum function.
#
# Arguments:
# beta_eta: a vector that concatenates the beta coefficient vector and 
# the eta coefficient vector, the point at which the log likelihood 
# is evaluated
# Y_name: a string denoting the outcome variable
# X_name: a string denoting the error-free exposure variable
# Z_name: a string denoting an additional error-free covariate 
# (optional, defaults to NULL)
# Q_name: a string denoting the query indicator variable
# offset_name: a string denoting the name of the offset variable 
# (optional, defaults to NULL)
# data: a data frame containing at least the columns specified above
# verbose: a logical parameter controlling the verbose functionality 
# (optional, defaults to FALSE)

loglik_mat = function(beta_eta, 
                      Y_name, X_name, 
                      Z_name = NULL, Xstar_name, 
                      Q_name, 
                      offset_name = NULL,
                      data,
                      verbose = FALSE) {
  
  # Save useful constants
  N = nrow(data) ## Phase I sample size
  n = sum(data[, Q_name]) ## Phase II sample size
  
  # Reorder data to put queried rows first
  data = data[order(data[, Q_name], decreasing = TRUE), ]
  
  # Create matrix of complete data
  if(!is.null(Z_name)){ #case with covariates
    if (n < N) {
      queried_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, 
                                                 Z_name, 
                                                 Xstar_name, 
                                                 offset_name)])
      unqueried_data = rbind(
        cbind(id = (n+1):N, data[-c(1:n), Y_name], 
              X_name = 0, data[-c(1:n), 
                               c(Z_name, Xstar_name, offset_name)]),
        cbind(id = (n+1):N, data[-c(1:n), Y_name], 
              X_name = 1, data[-c(1:n), 
                               c(Z_name, Xstar_name, offset_name)])
      )
      colnames(unqueried_data) = c("id", Y_name, X_name, Z_name, 
                                   Xstar_name, offset_name)
      complete_data = data.matrix(rbind(queried_data, unqueried_data))
    } else {
      complete_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, 
                                                  Z_name, 
                                                  Xstar_name, 
                                                  offset_name)])
    }
  } else{ #case without covariates
    if (n < N) {
      queried_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, 
                                                 Z_name, 
                                                 Xstar_name, 
                                                 offset_name)])
      unqueried_data = rbind(
        cbind(id = (n+1):N, data[-c(1:n), Y_name], 
              X_name = 0, data[-c(1:n), 
                               c(Xstar_name, offset_name)]),
        cbind(id = (n+1):N, data[-c(1:n), Y_name], 
              X_name = 1, data[-c(1:n), 
                               c(Xstar_name, offset_name)])
      )
      colnames(unqueried_data) = c("id", Y_name, X_name, 
                                   Xstar_name, offset_name)
      complete_data = data.matrix(rbind(queried_data, unqueried_data))
    } else {
      complete_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, 
                                                  Xstar_name, 
                                                  offset_name)])
    }
  }
  
  # Compute log-likelihood 
  if(!is.null(Z_name)){ ## P(Y|X,Z) from Poisson distribution
    lambdaY = exp(beta_eta[1] + beta_eta[2] * complete_data[, X_name] + 
                    beta_eta[3] * complete_data[, Z_name])  
    if(!is.null(offset_name)){ #has offset
      lambdaY = complete_data[, offset_name] * lambdaY
    } 
  } else{ ## P(Y|X) from Poisson distribution
    lambdaY = exp(beta_eta[1] + beta_eta[2] * complete_data[, X_name])
    if(!is.null(offset_name)){ #has offset
      lambdaY = complete_data[, offset_name] * lambdaY
    } 
  }

  pYgivXZ = dpois(x = complete_data[, Y_name], lambda = lambdaY)
  
  if(!is.null(Z_name)){ ## P(X|X*,Z) from Bernoulli distribution
    pXgivXstarZ = 1 / 
      (1 + 
         exp(-(beta_eta[4] + 
                 beta_eta[5] * complete_data[, Xstar_name] + 
                 beta_eta[6] * complete_data[, Z_name]))) ^ 
      complete_data[, X_name] * 
      (1 - 1 / (1 + exp(-(beta_eta[4] + 
                            beta_eta[5] * 
                            complete_data[, Xstar_name] + 
                            beta_eta[6] * 
                            complete_data[, Z_name])))) ^ 
      (1 - complete_data[, X_name]) 
  } else{ ## P(X|X*) from Bernoulli distribution
    pXgivXstarZ = 1 / 
      (1 + exp(-(beta_eta[3] + 
                   beta_eta[4] * complete_data[, Xstar_name]))) ^ 
      complete_data[, X_name] * 
      (1 - 1 / (1 + exp(-(beta_eta[3] +
                            beta_eta[4] * 
                            complete_data[, Xstar_name])))) ^ 
      (1 - complete_data[, X_name]) 
  }
  ## P(Y, X|X*, Z) OR P(Y,X|X*) 
  pYXgivXstarZ = pYgivXZ * pXgivXstarZ
  
  ## Marginalize X out of P(Y, X|X*, Z) for unqueried 
  marg_pYXgivXstarZ = rowsum(x = pYXgivXstarZ, 
                             group = complete_data[, "id"])
  
  ### replace zeroes with a VERY small number that's close to 0
  ### to avoid undefined logs
  pYgivXZ[which(pYgivXZ == 0)] = 5e-324
  pXgivXstarZ[which(pXgivXstarZ == 0)] = 5e-324
  marg_pYXgivXstarZ[which(marg_pYXgivXstarZ == 0)] = 5e-324
  
  # Compute log-likelihood
  ll = sum(log(pYgivXZ[c(1:n)])) + 
    sum(log(pXgivXstarZ[c(1:n)])) + 
    sum(log(marg_pYXgivXstarZ[-c(1:n)]))
  if(verbose) {print(paste("Queried:", ll))}
  
  ll = ll + 
    sum(log(marg_pYXgivXstarZ[-c(1:n)]))
  if(verbose) {print(paste("Queried + Unqueried:", ll))}
  return(-ll) ## return (-1) x log-likelihood for maximization
}

## Example
library(dplyr) #for data wrangling
set.seed(1031) #for reproducibility

#generate data
beta <- c(-2.2, 0.15) #governs Poisson outcome
eta <- c(-2.2, 4.4) #governs logistic error model
xstar = rbinom(n = 500, size = 1, prob = 0.5) #error-prone exposure
x = rbinom(n = 500, size = 1, #error-free exposure X|X*
             prob = 1 / (1 + exp(-(eta[1] + eta[2] * xstar)))) 
lambda = exp(beta[1] + beta[2] * x) #mean of Y|X
y = rpois(n = 500, lambda = lambda) #Poisson outcome with mean lambda
q = rbinom(n = 500, size = 1, prob = 0.75) ## queried indicator
df <- data.frame(xstar, x, y, q)
df <- df |> mutate(x = ifelse(q == 1, x, NA)) #redact X for unqueried rows

#call MLE function
mle_output <- mlePossum(error_formula = x ~ xstar,
                        analysis_formula = y ~ x,
                        data = df)
