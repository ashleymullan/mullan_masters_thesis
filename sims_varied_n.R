#Simulation Setup:
#Testing Our Estimator against Complete Case, Naive, and Gold Standard
#Varying Sample Size

library(dplyr)

eta_finder <- function(fpr, tpr){
  eta0 <- -1 * log((1 / fpr) - 1)
  eta1 <- -1 * (eta0 + log((1 / tpr) - 1)) 
  eta2 <- 0.5 #something I picked
  return(c(eta0, eta1, eta2))
}

loglik <- function(beta_eta, data, y, x, z, xstar, q, 
                   negate = TRUE, verbose = FALSE) {
  y <- data |> pull(y)
  x <- data |> pull(x)
  z <- data |> pull(z)
  xstar <- data |> pull(xstar)
  q <- data |> pull(q)
  ## If q = 1, return log(P(Y|X,Z)P(X|X*,Z))
  lambdaY = exp(beta_eta[1] + beta_eta[2] * x + beta_eta[3] * z)
  sigX = 1 / (1 + exp(- (beta_eta[4] + beta_eta[5] * xstar + beta_eta[6] * z)))
  inside_piece <- dpois(x = y, lambda = lambdaY) * sigX ^ x * (1 - sigX) ^ (1 - x)
  fixed <- ifelse(inside_piece == 0, 1, inside_piece)
  ll_q1 = sum(log(fixed), na.rm = TRUE)
  if(verbose) {print(paste("Queried contribution:", ll_q1))}
  
  ## If q = 0, return log(P(Y|X=0,Z)P(X=0|X*,Z) + P(Y|X=1,Z)P(X=1|X*,Z))
  lambdaX0 = exp(beta_eta[1] + beta_eta[2] * 0 + beta_eta[3] * z[q == 0]) ## vector of length (N-n)
  lambdaX1 = exp(beta_eta[1] + beta_eta[2] * 1 + beta_eta[3] * z[q == 0]) ## vector of length (N-n)
  sigX = 1 / (1 + exp(- (beta_eta[4] + beta_eta[5] * xstar[q == 0] + beta_eta[6] * z[q == 0]))) ## vector of length (N-n)
  inside_piece <- dpois(x = y[q == 0], lambda = lambdaX0) * sigX ^ 0 * (1 - sigX) ^ (1 - 0) + ## P(Y|X=0,Z)P(X=0|X*,Z)
    dpois(x = y[q == 0], lambda = lambdaX1) * sigX ^ 1 * (1 - sigX) ^ (1 - 1)
  fixed <- ifelse(inside_piece == 0, 1, inside_piece)
  ll_q0 = sum(log(fixed)) ## P(Y|X=1,Z)P(X=1|X*,Z)
  if(verbose) {print(paste("Unqueried contribution:", ll_q0))}
  ## Return sum of their contributions
  return((-1) ^ negate * (ll_q1 + ll_q0))
}

set.seed(1031) ## be a reproducible queen
num_sims <- 10000
results <- expand.grid(sim_id = 1:num_sims,
                       n = c(100, 1000, 10000), #medium, large, larger sample size
                       error_code = "M",
                       missingness = 0.50, 
                       betahat0_gs = NA, betahat1_gs = NA, betahat2_gs = NA, 
                       se_betahat0_gs = NA, se_betahat1_gs = NA, se_betahat2_gs = NA,
                       betahat0_n = NA, betahat1_n = NA, betahat2_n = NA, 
                       se_betahat0_n = NA, se_betahat1_n = NA, se_betahat2_n = NA,
                       betahat0_cc = NA, betahat1_cc = NA, betahat2_cc = NA, 
                       se_betahat0_cc = NA, se_betahat1_cc = NA, se_betahat2_cc = NA,
                       betahat0_my = NA, betahat1_my = NA, betahat2_my = NA,
                       se_betahat0_my = NA, se_betahat1_my = NA, se_betahat2_my = NA
)

results <- results |>
  mutate(fpr =
           case_when(error_code == "S" ~ 0.1, 
                     error_code == "M" ~ 0.25,
                     error_code == "L" ~ 0.5),
         .after = error_code
  ) |>
  mutate(tpr =
           case_when(error_code == "S" ~ 0.9, 
                     error_code == "M" ~ 0.75,
                     error_code == "L" ~ 0.5),
         .after = fpr
  ) 

for (r in 1:nrow(results)) { 
  # Simulate data 
  Z = rbinom(n = results$n[r], size = 1, prob = 0.3) ## Z ~ Bern(p = 0.3)
  Xstar = rbinom(n = results$n[r], size = 1, prob = 1 / (1 + exp(- (1 + 2 * Z))))
  e = eta_finder(results$fpr[r], results$tpr[r])
  beta = c(1,2,4)
  X = rbinom(n = results$n[r], size = 1, prob = 1 / (1 + exp(-(e[1] + e[2] * Xstar + e[3] * Z))))
  lambda = exp(beta[1] + beta[2] * X + beta[3] * Z) ## mean of the Poisson distribution for Y|X,Z
  Y = rpois(n = results$n[r], lambda = lambda) ## Y|X,Z ~ Pois(lambda), where lambda is a function of X, Z
  Q1 = rep(1, times = ceiling((1 - results$missingness[r]) * results$n[r]))
  Q0 = rep(0, times = results$n[r] - length(Q1))
  Q = c(Q1, Q0) ## queried indicator
  dat = data.frame(Y, X, Z, Xstar, Q)
  
  # Fit "gold standard" model
  gs = glm(formula = Y ~ X + Z, 
           family = poisson(link = "log"), 
           data = dat)
  
  # Save coefficient estimates from gold standard
  results[r, c("betahat0_gs", "betahat1_gs", "betahat2_gs")] = coefficients(gs)
  
  # Save standard error estimates from gold standard
  results[r, c("se_betahat0_gs", "se_betahat1_gs", "se_betahat2_gs")] = sqrt(diag(vcov(gs)))
  
  #fit naive model
  naive <- glm(Y ~ Xstar + Z, data = dat, family = poisson(link = "log"))
  
  # Save coefficient estimates from naive
  results[r, c("betahat0_n", "betahat1_n", "betahat2_n")] = coefficients(naive)
  
  # Save standard error estimates from naive
  results[r, c("se_betahat0_n", "se_betahat1_n", "se_betahat2_n")] = sqrt(diag(vcov(naive)))
  
  #fit complete case model
  complete <- dat |> filter(Q == 1)
  cc <- glm(Y ~ X + Z, data = complete, family = poisson(link = "log"))
  
  # Save coefficient estimates from complete case
  results[r, c("betahat0_cc", "betahat1_cc", "betahat2_cc")] = coefficients(cc)
  
  # Save standard error estimates from complete case
  results[r, c("se_betahat0_cc", "se_betahat1_cc", "se_betahat2_cc")] = sqrt(diag(vcov(cc)))
  
  #fit our estimator
  ours <- optim(par = rep(0, times = 6),
                fn = loglik,
                # gr = ??, #FILL ME IN
                data = dat,
                y = "Y",
                x = "X",
                z = "Z",
                xstar = "Xstar",
                q = "Q",
                method = "BFGS",
                hessian = TRUE)
  
  # Save coefficient estimates from our estimator
  results[r, c("betahat0_my", "betahat1_my", "betahat2_my")] = ours$par[1:3]
  
  # Save standard error estimates from our estimator
  #results[r, c("se_betahat0_my", "se_betahat1_my", "se_betahat2_my")] = sqrt(diag(solve(ours$hessian)))
  if(r %% 100 == 0) {print(paste("run number ", r))}
}

results |> group_by(n) |> summarize(
  betahat0_cc = mean(betahat0_cc, na.rm = TRUE),
  betahat1_cc = mean(betahat1_cc, na.rm = TRUE),
  betahat2_cc = mean(betahat2_cc, na.rm = TRUE),
  betahat0_gs = mean(betahat0_gs, na.rm = TRUE),
  betahat1_gs = mean(betahat1_gs, na.rm = TRUE),
  betahat2_gs = mean(betahat2_gs, na.rm = TRUE),
  betahat0_n = mean(betahat0_n, na.rm = TRUE),
  betahat1_n = mean(betahat1_n, na.rm = TRUE),
  betahat2_n = mean(betahat2_n, na.rm = TRUE),
  betahat0_my = mean(betahat0_my, na.rm = TRUE),
  betahat1_my = mean(betahat1_my, na.rm = TRUE),
  betahat2_my = mean(betahat2_my, na.rm = TRUE),
) |> View()
