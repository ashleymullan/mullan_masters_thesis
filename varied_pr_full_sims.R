eta_finder <- function(fpr, tpr){
  # Parameters for error model P(X*|X,Z) ----------------------
  eta0 <- - log((1 - fpr) / fpr) 
  eta1 <- - log((1 - tpr) / tpr) - eta0
  eta2 <- 0.5
  return(c(eta0, eta1, eta2))
}

# True values of the model coefficients 
b0 = 2 ## intercept
b1 = 3 ## log prevalence ratio for X (conditioning on Z)
b2 = 4 ## log prevalence ratio for Z (conditioning on X)
e = eta_finder(fpr = 0.1, tpr = 0.9)

loglik_mat = function(beta_eta, 
                      Y_name, X_name, 
                      Z_name, Xstar_name, 
                      Q_name, data,
                      verbose = FALSE) {
  #print(beta_eta)
  
  # Save useful constants
  N = nrow(data) ## Phase I sample size
  n = sum(data[, Q_name]) ## Phase II sample size
  
  # Reorder data to put queried rows first
  data = data[order(data[, Q_name], decreasing = TRUE), ]
  
  # Create matrix of complete data
  if (n < N) {
    queried_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name)])
    unqueried_data = rbind(
      cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 0, data[-c(1:n), c(Z_name, Xstar_name)]),
      cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 1, data[-c(1:n), c(Z_name, Xstar_name)])
    )
    colnames(unqueried_data) = c("id", Y_name, X_name, Z_name, Xstar_name)
    complete_data = data.matrix(rbind(queried_data, unqueried_data))
  } else {
    complete_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name)])
  }
  
  # Compute log-likelihood 
  ## P(Y|X,Z) from Poisson distribution
  lambdaY = exp(beta_eta[1] + beta_eta[2] * complete_data[, X_name] + beta_eta[3] * complete_data[, Z_name])
  
  ### Dazzle fix: replace y with data[, Y_name]
  pYgivXZ = dpois(x = complete_data[, Y_name], lambda = lambdaY)
  
  ## P(X|X*,Z) from Bernoulli distribution
  pXgivXstarZ = 1 / (1 + exp(-(beta_eta[4] + beta_eta[5] * complete_data[, Xstar_name] + beta_eta[6] * complete_data[, Z_name]))) ^ complete_data[, X_name] * 
    (1 - 1 / (1 + exp(-(beta_eta[4] + beta_eta[5] * complete_data[, Xstar_name] + beta_eta[6] * complete_data[, Z_name])))) ^ (1 - complete_data[, X_name]) 
  
  ## P(Y, X|X*, Z) 
  pYXgivXstarZ = pYgivXZ * pXgivXstarZ
  
  ## Marginalize X out of P(Y, X|X*, Z) for unqueried 
  marg_pYXgivXstarZ = rowsum(x = pYXgivXstarZ, 
                             group = complete_data[, "id"])
  
  ### Dazzle fix: replace with another VERY small number that's close to 0
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

# Simulation to check that the "gold standard" model returns correct estimates 
set.seed(1031) ## be a reproducible queen
num_reps = 18000

res = data.frame(rep = 1:num_reps, code = NA, 
                 our_beta0 = NA, our_beta1 = NA, our_beta2 = NA, 
                 cc_beta0 = NA, cc_beta1 = NA, cc_beta2 = NA,
                 naive_beta0 = NA, naive_beta1 = NA, naive_beta2 = NA,
                 gs_beta0 = NA, gs_beta1 = NA, gs_beta2 = NA,
                 eta0 = NA, eta1 = NA, eta2 = NA,
                 n = c(rep(100, times = num_reps / 3), 
                       rep(1000, times = num_reps / 3), 
                       rep(10000, times = num_reps / 3)),
                 q = rep(0.75, times = num_reps),
                 fpr = rep(c(rep(0.1, times = num_reps / 9), 
                       rep(0.25, times = num_reps / 9), 
                       rep(0.5, times = num_reps / 9)), times = 3),
                 tpr = rep(c(rep(0.9, times = num_reps / 9), 
                       rep(0.75, times = num_reps / 9), 
                       rep(0.5, times = num_reps / 9)), times = 3),
                 our_beta0_se = NA, our_beta1_se = NA, our_beta2_se = NA, 
                 cc_beta0_se = NA, cc_beta1_se = NA, cc_beta2_se = NA,
                 naive_beta0_se = NA, naive_beta1_se = NA, naive_beta2_se = NA,
                 gs_beta0_se = NA, gs_beta1_se = NA, gs_beta2_se = NA)
print(paste("current time:", Sys.time()))
for (r in 1:num_reps) {
  # Simulate data 
  e = eta_finder(fpr = res$fpr[r], tpr = res$tpr[r])
  z = rnorm(n = res$n[r]) #rbinom(n = 10000, size = 1, prob = 0.3) ## Z ~ Bern(p = 0.3)
  xstar = rbinom(n = res$n[r], size = 1, prob = 1 / (1 + exp(- (1 + 2 * z)))) 
  x = rbinom(n = res$n[r], size = 1, prob = 1 / (1 + exp(-(e[1] + e[2] * xstar + e[3] * z))))
  lambda = exp(b0 + b1 * x + b2 * z) ## mean of the Poisson distribution for Y|X,Z
  y = rpois(n = res$n[r], lambda = lambda) ## Y|X,Z ~ Pois(lambda), where lambda is a function of X, Z
  q = rbinom(n = res$n[r], size = 1, prob = res$q[r]) #0.25) ## queried indicator
  dat = data.frame(y, x, z, xstar, q) 
  
  cc = glm(formula = y ~ x + z, 
           data = dat, 
           family = poisson, 
           subset = q == 1)
  
  cc_se = summary(cc)$coefficients[,"Std. Error"]
  
  cc_fit = glm(formula = y ~ x + z, 
               data = dat, 
               family = poisson, 
               subset = q == 1)$coefficients
  
  cc_fit = c(cc_fit, 
             glm(formula = x ~ xstar + z, 
                 data = dat, 
                 family = binomial, 
                 subset = q == 1)$coefficients)
  
  naive_fit = summary(glm(formula = y ~ xstar + z,
                          data = dat,
                          family = poisson))$coefficients
  
  gs_fit = summary(glm(formula = y ~ x + z,
                       data = dat,
                       family = poisson))$coefficients
  
  optim_res = optim(fn = loglik_mat, 
                    par = cc_fit, 
                    hessian = TRUE, 
                    method = "BFGS",
                    Y_name = "y",
                    X_name = "x",
                    Z_name = "z",
                    Xstar_name = "xstar",
                    Q_name = "q",
                    data = dat)
  num_analysis_covar = 3
  optim_vcov = tryCatch(expr = solve(optim_res$hessian)[1:num_analysis_covar, 1:num_analysis_covar],
                        error = function(err) {
                          matrix(data = NA, 
                                 nrow = num_analysis_covar, 
                                 ncol = num_analysis_covar)
                        })
  
  res[r, "code"] = optim_res$convergence
  res[r, "our_beta0"] = optim_res$par[1]
  res[r, "our_beta1"] = optim_res$par[2]
  res[r, "our_beta2"] = optim_res$par[3]
  res[r, "our_beta0_se"] = sqrt(diag(optim_vcov))[1]
  res[r, "our_beta1_se"] = sqrt(diag(optim_vcov))[2]
  res[r, "our_beta2_se"] = sqrt(diag(optim_vcov))[3]
  res[r, "eta0"] = optim_res$par[4]
  res[r, "eta1"] = optim_res$par[5]
  res[r, "eta2"] = optim_res$par[6]
  res[r, "cc_beta0"] = cc_fit[1]
  res[r, "cc_beta1"] = cc_fit[2]
  res[r, "cc_beta2"] = cc_fit[3]
  res[r, "cc_beta0_se"] = cc_se[1]
  res[r, "cc_beta1_se"] = cc_se[2]
  res[r, "cc_beta2_se"] = cc_se[3]
  res[r, "naive_beta0"] = naive_fit[1,"Estimate"]
  res[r, "naive_beta1"] = naive_fit[2,"Estimate"]
  res[r, "naive_beta2"] = naive_fit[3,"Estimate"]
  res[r, "naive_beta0_se"] = naive_fit[1,"Std. Error"]
  res[r, "naive_beta1_se"] = naive_fit[2,"Std. Error"]
  res[r, "naive_beta2_se"] = naive_fit[3,"Std. Error"]
  res[r, "gs_beta0"] = gs_fit[1, "Estimate"]
  res[r, "gs_beta1"] = gs_fit[2, "Estimate"]
  res[r, "gs_beta2"] = gs_fit[3, "Estimate"]
  res[r, "gs_beta0_se"] = gs_fit[1, "Std. Error"]
  res[r, "gs_beta1_se"] = gs_fit[2, "Std. Error"]
  res[r, "gs_beta2_se"] = gs_fit[3, "Std. Error"]
  
  if(r %% 100 == 0) {
    print(paste("rep #:", r))
    print(paste("current time:", Sys.time()))}
}

write.csv(res, "varied_pr_full_sims.csv")

