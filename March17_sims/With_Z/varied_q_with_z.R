library(devtools)
devtools::install_github(repo = "sarahlotspeich/possum")
library(possum)

eta_finder <- function(fpr, tpr){
  # Parameters for error model P(X*|X) ----------------------
  eta0 <- - log((1 - fpr) / fpr) 
  eta1 <- - log((1 - tpr) / tpr) - eta0
  eta2 <- 0.5
  return(c(eta0, eta1, eta2))
}

# True values of the model coefficients 
b0 = -2.246 ## intercept, real data gold standard fit
b1 = 0.155 ## log prevalence ratio for X, real data gold standard fit
b2 = 0.2
e = eta_finder(fpr = 0.1, tpr = 0.9) 


set.seed(1031) ## be a reproducible queen
num_reps = 16000
res = data.frame(rep = 1:num_reps, 
                 our_beta0 = NA, our_beta1 = NA, our_beta2 = NA,
                 cc_beta0 = NA, cc_beta1 = NA, cc_beta2 = NA,
                 naive_beta0 = NA, naive_beta1 = NA, naive_beta2 = NA,
                 gs_beta0 = NA, gs_beta1 = NA, gs_beta2 = NA,
                 n = c(rep(390, times = num_reps / 2), 
                       rep(2200, times = num_reps / 2)),
                 q = rep(c(rep(0.25, times = num_reps / 8),
                           rep(0.50, times = num_reps / 8), 
                           rep(0.75, times = num_reps / 8), 
                           rep(0.90, times = num_reps / 8)), times = 2),
                 our_beta0_se = NA, our_beta1_se = NA, our_beta2_se = NA,
                 cc_beta0_se = NA, cc_beta1_se = NA, cc_beta2_se = NA,
                 naive_beta0_se = NA, naive_beta1_se = NA, naive_beta2_se = NA,
                 gs_beta0_se = NA, gs_beta1_se = NA, gs_beta2_se = NA)


for (r in 1:num_reps) {
  # Simulate data 
  z = rnorm(n = res$n[r])
  xstar = rbinom(n = res$n[r], size = 1, prob = 1 / (1 + exp(- (1 + 2 * z)))) #chosen by me
  x = rbinom(n = res$n[r], size = 1, prob = 1 / (1 + exp(-(e[1] + e[2] * xstar + e[3] * z))))
  lambda = exp(b0 + b1 * x + b2 * z) ## mean of the Poisson distribution for Y|X,Z
  y = rpois(n = res$n[r], lambda = lambda) ## Y|X,Z ~ Pois(lambda), where lambda is a function of X
  q = rbinom(n = res$n[r], size = 1, prob = res$q[r]) ## queried indicator
  dat = data.frame(y, x, z, xstar, q) 
  
  cc = summary(glm(formula = y ~ x + z, 
                   data = dat, 
                   family = poisson, 
                   subset = q == 1))$coefficients
  
  naive = summary(glm(formula = y ~ xstar + z,
                      data = dat,
                      family = poisson))$coefficients
  
  gs = summary(glm(formula = y ~ x + z,
                   data = dat,
                   family = poisson))$coefficients
  
  dat$x <- ifelse(dat$q == 1, dat$x, NA) #set up for use in mlePossum
  
  mle = mlePossum(error_formula = x ~ xstar + z,
                  analysis_formula = y ~ x + z,
                  data = dat)

  res[r, "our_beta0"] = mle[1,1]
  res[r, "our_beta1"] = mle[2,1]
  res[r, "our_beta2"] = mle[3,1]
  res[r, "our_beta0_se"] = mle[1,2]
  res[r, "our_beta1_se"] = mle[2,2]
  res[r, "our_beta2_se"] = mle[3,2]
  res[r, "cc_beta0"] = cc[1,"Estimate"]
  res[r, "cc_beta1"] = cc[2,"Estimate"]
  res[r, "cc_beta2"] = cc[3,"Estimate"]
  res[r, "cc_beta0_se"] = cc[1,"Std. Error"]
  res[r, "cc_beta1_se"] = cc[2,"Std. Error"]
  res[r, "cc_beta2_se"] = cc[3,"Std. Error"]
  res[r, "naive_beta0"] = naive[1,"Estimate"]
  res[r, "naive_beta1"] = naive[2,"Estimate"]
  res[r, "naive_beta2"] = naive[3,"Estimate"]
  res[r, "naive_beta0_se"] = naive[1,"Std. Error"]
  res[r, "naive_beta1_se"] = naive[2,"Std. Error"]
  res[r, "naive_beta2_se"] = naive[3,"Std. Error"]
  res[r, "gs_beta0"] = gs[1, "Estimate"]
  res[r, "gs_beta1"] = gs[2, "Estimate"]
  res[r, "gs_beta2"] = gs[3, "Estimate"]
  res[r, "gs_beta0_se"] = gs[1, "Std. Error"]
  res[r, "gs_beta1_se"] = gs[2, "Std. Error"]
  res[r, "gs_beta2_se"] = gs[3, "Std. Error"]

  
}

write.csv(res, "vqwz.csv")

