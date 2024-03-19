library(devtools)
devtools::install_github(repo = "sarahlotspeich/possum")
library(possum)

eta_finder <- function(fpr, tpr){
  # Parameters for error model P(X*|X) ----------------------
  eta0 <- - log((1 - fpr) / fpr) 
  eta1 <- - log((1 - tpr) / tpr) - eta0
  return(c(eta0, eta1))
}

# True values of the model coefficients 
b0 = -2.246
e = eta_finder(fpr = 0.1, tpr = 0.9) 


set.seed(1031) ## be a reproducible queen
num_reps = 24000
res = data.frame(rep = 1:num_reps, 
                 our_beta0 = NA, our_beta1 = NA, 
                 cc_beta0 = NA, cc_beta1 = NA, 
                 naive_beta0 = NA, naive_beta1 = NA,
                 gs_beta0 = NA, gs_beta1 = NA,
                 q = rep(0.75, times = num_reps),
                 n = c(rep(390, times = num_reps / 2), 
                       rep(2200, times = num_reps / 2)),
                 b1 = rep(c(rep(-5, times = num_reps / 8),
                           rep(0.155, times = num_reps / 8), 
                           rep(10, times = num_reps / 8),
                           rep(20, times = num_reps / 8)), times = 2),
                 our_beta0_se = NA, our_beta1_se = NA,  
                 cc_beta0_se = NA, cc_beta1_se = NA, 
                 naive_beta0_se = NA, naive_beta1_se = NA, 
                 gs_beta0_se = NA, gs_beta1_se = NA)

for (r in 1:num_reps) {
  # Simulate data 
  b1 = res$b1[r]
  xstar = rbinom(n = res$n[r], size = 1, prob = 0.496) #from real data prevalence at 1 mile level 
  x = rbinom(n = res$n[r], size = 1, prob = 1 / (1 + exp(-(e[1] + e[2] * xstar))))
  lambda = exp(b0 + b1 * x) ## mean of the Poisson distribution for Y|X
  y = rpois(n = res$n[r], lambda = lambda) ## Y|X ~ Pois(lambda), where lambda is a function of X
  q = rbinom(n = res$n[r], size = 1, prob = res$q[r]) ## queried indicator
  dat = data.frame(y, x, xstar, q) 
  
  cc = summary(glm(formula = y ~ x, 
                   data = dat, 
                   family = poisson, 
                   subset = q == 1))$coefficients
  
  naive = summary(glm(formula = y ~ xstar,
                      data = dat,
                      family = poisson))$coefficients
  
  gs = summary(glm(formula = y ~ x,
                   data = dat,
                   family = poisson))$coefficients
  
  dat$x <- ifelse(dat$q == 1, dat$x, NA) #set up for use in mlePossum
  
  mle = mlePossum(error_formula = x ~ xstar,
                  analysis_formula = y ~ x,
                  data = dat)
  
  res[r, "our_beta0"] = mle[1,1]
  res[r, "our_beta1"] = mle[2,1]
  res[r, "our_beta0_se"] = mle[1,2]
  res[r, "our_beta1_se"] = mle[2,2]
  res[r, "cc_beta0"] = cc[1,"Estimate"]
  res[r, "cc_beta1"] = cc[2,"Estimate"]
  res[r, "cc_beta0_se"] = cc[1,"Std. Error"]
  res[r, "cc_beta1_se"] = cc[2,"Std. Error"]
  res[r, "naive_beta0"] = naive[1,"Estimate"]
  res[r, "naive_beta1"] = naive[2,"Estimate"]
  res[r, "naive_beta0_se"] = naive[1,"Std. Error"]
  res[r, "naive_beta1_se"] = naive[2,"Std. Error"]
  res[r, "gs_beta0"] = gs[1, "Estimate"]
  res[r, "gs_beta1"] = gs[2, "Estimate"]
  res[r, "gs_beta0_se"] = gs[1, "Std. Error"]
  res[r, "gs_beta1_se"] = gs[2, "Std. Error"]
  
}

write.csv(res, "vb1nz.csv")

