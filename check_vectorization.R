eta_finder <- function(fpr, tpr){
  eta0 <- -1 * log((1 / fpr) - 1)
  eta1 <- -1 * (eta0 + log((1 / tpr) - 1)) 
  eta2 <- 0.5 #something I picked
  return(c(eta0, eta1, eta2))
}

generate_data <- function(n, a, fpr, tpr, beta, qp){
  eta <- eta_finder(fpr, tpr)
  Z <- rbinom(n, 1, 0.2) #fake rural indicator
  
  z_dependence <- sf(a[1] + a[2]*Z) # P(X* = 1|Z)
  Xstar <- rbinom(n, 1, z_dependence) # error-prone food access indicator
  
  zxstar_dependence <- sf(eta[1] + eta[2]*Xstar + eta[3]*Z) # P(X=1|X*,Z)
  X <- rbinom(n, 1, zxstar_dependence)
  
  lambda <- exp(beta[1] + beta[2]*X + beta[3]*Z)
  Y <- rpois(n, lambda)
  Q1 <- rep(1, times = ceiling(qp * n)) #picks the first query percentage elements to be queried
  #note to user: this is fine because we already generated data randomly!
  Q0 <- rep(0, times = n - length(Q1))
  Q <- c(Q1, Q0)
  return(data.frame(X, Xstar, Y, Z, Q))
}

#sigmoid function
sf <- function(x){
  1 / (1 + exp(-x))
}


negative_ell <- function(beta_eta, data){
  
  beta <- beta_eta[1:3]
  eta <- beta_eta[4:6]
  
  X <- data$X #length N
  Xstar <- data$Xstar 
  Y <- data$Y
  Z <- data$Z
  Q <- data$Q
  
  N <- nrow(data)
  n <- sum(Q)
  
  queried <- data |> filter(Q == 1)
  q_contribution <- 0 #initialize queried contribution
  unqueried <- data |> filter(Q == 0)
  uq_contribution <- 0 #initialize unqueried contribution
  
  #alternative option for contributing to queried likelihood
  for(i in 1:nrow(queried)){
    y_given_xz <- dpois(queried$Y[i], exp(beta[1] + beta[2]*queried$X[i] + beta[3]*queried$Z[i]))
    p <- sf(eta[1] + eta[2]*queried$Xstar[i] + eta[3]*queried$Z[i]) #p in binomial, computed as logistic
    x_given_xstarz <- p^queried$X[i] * (1-p)^(1 - queried$X[i])
    #p(x*,z) is proportional and drops out
    lik <- y_given_xz * x_given_xstarz
    q_contribution <- q_contribution + ifelse(lik != 0, log(lik), 0)
  }
  
  #alternative option for contributing to unqueried likelihood
  for(i in 1:nrow(unqueried)){
    #do separately for both Xs, similar style to the loop starting on 8
    #X = 0
    y_given_xz <- dpois(unqueried$Y[i], exp(beta[1] + beta[3]*unqueried$Z[i]))
    p <- sf(eta[1] + eta[2]*unqueried$Xstar[i] + eta[3]*unqueried$Z[i])
    x_given_xstarz <- (1- p) #X is zero so I didn't bother writing the other piece of the product
    #p(x*,z) is proportional and drops out
    lik_x0 <- y_given_xz * x_given_xstarz
    
    #X = 1
    y_given_xz <- dpois(unqueried$Y[i], exp(beta[1] + beta[2] + beta[3]*unqueried$Z[i]))
    p <- sf(eta[1] + eta[2]*unqueried$Xstar[i] + eta[3]*unqueried$Z[i])
    x_given_xstarz <- p #X is 1 so I didn't bother writing the other piece of the product
    #p(x*,z) is proportional and drops out
    lik_x1 <- y_given_xz * x_given_xstarz
    lik <- lik_x0 + lik_x1
    uq_contribution <- uq_contribution + ifelse(lik != 0, log(lik), 0)
  }
  return(-1 *(q_contribution + uq_contribution))
}


num_sims <- 500 #FIX ME
a <- c(1,2) #parameters for Z dependence, generating X*|Z
beta <- c(1,2,3) #parameters of Y|X,Z

results <- expand.grid(sim_id = 1:num_sims,
                       n = c(100, 1000, 10000), #medium, large, larger sample size
                       error_code = c("S", "M", "L"),
                       missingness = c(0, 0.25, 0.50, 0.75), #naive, three levels of complete
                       betahat0_gs = NA, betahat1_gs = NA, betahat2_gs = NA, 
                       se_betahat0_gs = NA, se_betahat1_gs = NA, se_betahat2_gs = NA,
                       betahat0_n = NA, betahat1_n = NA, betahat2_n = NA, 
                       se_betahat0_n = NA, se_betahat1_n = NA, se_betahat2_n = NA,
                       betahat0_cc = NA, betahat1_cc = NA, betahat2_cc = NA, 
                       se_betahat0_cc = NA, se_betahat1_cc = NA, se_betahat2_cc = NA,
                       betahat0_my = NA, betahat1_my = NA, betahat2_my = NA
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

#Sims of Bad Methods and My Badly Written Method
set.seed(1031)
for (i in 1:nrow(results)) {
  df <- generate_data(n = results$n[i], 
                      a = a,
                      beta = beta,
                      fpr = results$fpr[i],
                      tpr = results$tpr[i],
                      qp = 0.5) #initialize simulated data
  
  eta <- eta_finder(df$fpr[i], df$tpr[i])
  
  
  #gold standard
  gs <- glm(Y ~ X + Z, data = df, family = "poisson")
  gs <- gs |> tidy()
  results[i, c("betahat0_gs", "betahat1_gs", "betahat2_gs")] <- gs |> pull(estimate)
  results[i, c("se_betahat0_gs", "se_betahat1_gs", "se_betahat2_gs")] <- gs |> pull(std.error)
  
  #naive
  naive <- glm(Y ~ Xstar + Z, data = df, family = "poisson")
  naive <- naive |> tidy()
  results[i, c("betahat0_n", "betahat1_n", "betahat2_n")] <- naive |> pull(estimate)
  results[i, c("se_betahat0_n", "se_betahat1_n", "se_betahat2_n")] <- naive |> pull(std.error)
  
  #complete case
  complete <- df |> filter(Q == 1)
  cc <- glm(Y ~ X + Z, data = complete, family = "poisson")
  cc <- cc |> tidy()
  results[i, c("betahat0_cc", "betahat1_cc", "betahat2_cc")] <- cc |> pull(estimate)
  results[i, c("se_betahat0_cc", "se_betahat1_cc", "se_betahat2_cc")] <- cc |> pull(std.error)
  
  #my estimator
  ests <- nlm(f = negative_ell, 
              p = c(0,0,0,0,0,0),
              data = df)
  results[i, c("betahat0_my", "betahat1_my", "betahat2_my")] <- ests$estimate[1:3]
  
  #sim counter
  if (i %% 10 == 0) {
    print(i)
  }
}

results |> group_by(n, missingness, error_code) |>
  summarize(bh0 = mean(betahat0_my),
            bh1 = mean(betahat1_my),
            bh2 = mean(betahat2_my))

