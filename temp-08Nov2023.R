# Run once: 
## devtools::install_github(repo = "sarahlotspeich/possum")

set.seed(1031)

beta0 <- 0.5
beta1 <- 0.25
beta2 <- 0.2
x <- runif(100, min = 0, max = 1) #generate 100 rows from U(0,1)
z <- runif(100, min = 0.1, max = 1.1) #generate 100 rows from U(0.1,1.1)
lambda <- exp(beta0 + beta1*x + beta2*z)
y <- rpois(100, lambda)

possum::mlePossum(start_guess = c(0.1, 0.2, 0.21), 
                  x = "x", 
                  y = "y", 
                  z = "z", 
                  data = data.frame(x,y,z))

my_est <- round(find_beta(start_guess = c(0.1, 0.2, 0.21), 
                          x = "x", 
                          y = "y", 
                          z = "z", 
                          data = data.frame(x,y,z))$estimates,2) 

default_est <- round(glm(y ~ x + z, family = "poisson")$coefficients,2)