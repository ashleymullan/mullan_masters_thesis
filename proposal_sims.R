
library(dplyr)
library(ggplot2)
library(latex2exp)
library(broom)
library(ggthemes)
library(tidyr)

cbbPalette <- c("#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00", "#56B4E9", "#CC79A7", "#000000")


eta_finder <- function(fpr, tpr){
  eta0 <- -1 * log((1 / fpr) - 1)
  eta1 <- -1 * (eta0 + log((1 / tpr) - 1)) 
  eta2 <- 0.5 #something I picked
  return(c(eta0, eta1, eta2))
}

generate_data <- function(n, a, fpr, tpr, beta){
  eta <- eta_finder(fpr, tpr)
  Z <- rbinom(n, 1, 0.2) #fake rural indicator
  
  z_dependence <- 1 / (1 + exp(- (a[1] + a[2]*Z))) # P(X* = 1|Z)
  Xstar <- rbinom(n, 1, z_dependence) # error-prone food access indicator
  
  zxstar_dependence <- 1 / (1 + exp(- (eta[1] + eta[2]*Xstar + eta[3]*Z))) # P(X=1|X*,Z)
  X <- rbinom(n, 1, zxstar_dependence)
  
  lambda <- exp(beta[1] + beta[2]*X + beta[3]*Z)
  Y <- rpois(n, lambda)
  

  return(data.frame(X, Xstar,Y,Z))
}


#TASK: it runs, figure out why my coefficient estimates are huge
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
  
 
  #dealing with the queried contribution to the likelihood
  p_ygiven_xz <- dpois(Y[Q == 1], exp(beta[1] + beta[2]*X[Q == 1] + beta[3]*Z[Q == 1]))
  p1_Q1 <- 1 / (1 + exp(- (eta[1] + eta[2]*Xstar[Q == 1] + eta[3]*Z[Q == 1]))) #length n
  p_x_given_xstarz <- p1_Q1^X[Q == 1] * (1 - p1_Q1)^(1-X[Q == 1])
  
  P4 <- p_ygiven_xz * p_x_given_xstarz #length n
  
  #dealing with the non-queried contribution to the likelihood
  
  x_options <- matrix(data = c(0,1), nrow = N - n, ncol = 2, byrow = TRUE)
  
  P3 <- rep(0, times = N - n)
  
  for(i in 1:ncol(x_options)){
    
    p_ygiven_xz <- dpois(Y[Q == 0], exp(beta[1] + beta[2]*x_options[,i] + beta[3]*Z[Q == 0]))
    p1_Q1 <- 1 / (1 + exp(- (eta[1] + eta[2]*Xstar[Q == 0] + eta[3]*Z[Q == 0]))) #length n
    p_x_given_xtarz <- p1_Q1^x_options[,i] * (1 - p1_Q1)^(1-x_options[,i])  
    P3 <- P3 + p_ygiven_xz * p_x_given_xstarz

  }
    
    P3[P3 == 0] <- 1 #to stop the issue at log(0), yields zero contribution anyway
    P4[P4 == 0] <- 1 #to stop the issue at log(0), yields zero contribution anyway
    
  return(-1*(sum(log(P4)) + sum(log(P3))))
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
                       se_betahat0_cc = NA, se_betahat1_cc = NA, se_betahat2_cc = NA
                       
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

#Sims of Bad Methods
set.seed(1031)
for (i in 1:nrow(results)) {
  df <- generate_data(n = results$n[i], 
                      a = a,
                      beta = beta,
                      fpr = results$fpr[i],
                      tpr = results$tpr[i]) #initialize simulated data
  
  eta <- eta_finder(df$fpr[i], df$tpr[i])
  
  df <- df |>
    mutate(id = 1:nrow(df)) #assign row ID to simulated data
  
  kept_ids <- sample(df$id, 
                     size = (1 - results$missingness[i]) * nrow(df),
                     replace = FALSE)
  
  df <- df |>
    mutate(Q = ifelse(id %in% kept_ids, 1, 0)) #figure out Q
  
  
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
}

results_long <- results %>%
  pivot_longer(cols = starts_with("betahat1_"), 
               names_to = "method_type", 
               values_to = "betahat1")


#plot empirical distributions
fixed_qual_plot <- 
  results_long |>
  filter(missingness == 0.5) |> #fix missingness to 50%
  filter(error_code == "M") |> #fix medium error
  mutate(method_type = case_when(
    method_type == "betahat1_gs" ~ "Gold Standard", 
    method_type == "betahat1_cc" ~ "Complete Case",
    method_type == "betahat1_n" ~ "Naive")) |>
  #filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(n), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('100','1,000','10,000')) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 2, color = "gray", linetype = "dashed") +
  scale_fill_manual(values = cbbPalette, name = "Method") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       fill = "Missingness",
       x = "(Original) Sample Size",
       title = "Method Comparison with Fixed Data Quality",
       subtitle = "50% Missingness and Medium Error Setting")
ggsave(filename = "fixed_quality_results.pdf", 
       plot = fixed_qual_plot,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")


#plot empirical distributions
varying_error_plot <- 
  results_long |>
  filter(missingness == 0.5) |> #fix missingness to 50%
  filter(n == 10000) |> #fix medium sample size
  mutate(method_type = case_when(
    method_type == "betahat1_gs" ~ "Gold Standard", 
    method_type == "betahat1_cc" ~ "Complete Case",
    method_type == "betahat1_n" ~ "Naive")) |>
  filter(method_type != "Complete Case") |> #temporarily axe complete case
  #filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(error_code), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('Small','Medium','Large')) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 2, color = "gray", linetype = "dashed") +
  scale_fill_manual(values = cbbPalette, name = "Method") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       fill = "Missingness",
       x = "Level of Error",
       title = "Method Comparison with Fixed Sample Size",
       subtitle = "N = 10,000")
ggsave(filename = "varying_error_results.png", 
       plot = varying_error_plot,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")

#plot empirical distributions
varying_missingness_plot <- 
  results_long |>
  filter(error_code == "M") |> #fix medium error
  filter(n == 10000) |> #fix large sample size
  filter(missingness > 0) |>
  mutate(method_type = case_when(
    method_type == "betahat1_gs" ~ "Gold Standard", 
    method_type == "betahat1_cc" ~ "Complete Case",
    method_type == "betahat1_n" ~ "Naive")) |>
  #filter(method_type != "Complete Case") |> #temporarily axe complete case
  filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(missingness), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('25%','50%', "75%")) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 2, color = "gray", linetype = "dashed") +
  scale_fill_manual(values = c(cbbPalette[2], cbbPalette[1]), name = "Method") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       x = "Missingness",
       title = "Method Comparison with Fixed Sample Size",
       subtitle = "N = 10,000")
ggsave(filename = "varying_missingness_results.png", 
       plot = varying_missingness_plot,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")


























#POST PROPOSAL SIMS

num_sims <- 500 #FIX ME
a <- c(1,2) #parameters for Z dependence, generating X*|Z
beta <- c(1,2,3) #parameters of Y|X,Z

results <- expand.grid(sim_id = 1:num_sims,
                       n = c(100, 1000, 10000), #medium, large, larger sample size
                       fpr = c(0.1, 0.25, 0.5), #small, medium, large error
                       tpr = c(0.9, 0.75, 0.5), #small, medium, large error
                       missingness = c(0, 0.25, 0.50, 0.75), #naive, three levels of complete
                       betahat0 = NA, betahat1 = NA, betahat2 = NA, 
                       converged = NA,
                       se_betahat0 = NA, se_betahat1 = NA, se_betahat2 = NA
)



# Set seed before looping
set.seed(1031)
for (i in 1:nrow(results)) {
  df <- generate_data(n = results$n[i], 
                      a = a,
                      beta = beta,
                      fpr = results$fpr[i],
                      tpr = results$tpr[i]) #initialize simulated data
  
  eta <- eta_finder(df$fpr[i], df$tpr[i])
  
  df <- df |>
    mutate(id = 1:nrow(df)) #assign row ID to simulated data
  
  kept_ids <- sample(df$id, 
                     size = (1 - results$missingness[i]) * nrow(df),
                     replace = FALSE)
  
  df <- df |>
    mutate(Q = ifelse(id %in% kept_ids, 1, 0)) #figure out Q
  
  
  optimized <- nlm(negative_ell, rep(0, times = 6), hessian = TRUE, data = df)
  estimates <- optimized$estimate
  convergence_code <- optimized$code
  #se <- sqrt(diag(solve(-optimized$hessian)))
  
  results[i, c("betahat0", "betahat1", "betahat2")] <- estimates
  results[i, c("converged")] <- convergence_code
  #results[i, c("se_betahat0", "se_betahat1", "se_betahat2")] <- se
}



results |>
  group_by(n, fpr, tpr, missingness) |>
  summarize(avg_coef = mean(betahat1),
            avg_se = mean(se_betahat1))

