
gradient <- function(dat, #data frame
                     xname, xstarname, zname, yname, qname, #strings
                     beta0, beta1, beta2, nu0, nu1, nu2){ #parameter guesses
  queried <- dat |> subset(dat[[qname]] == 1)
  unqueried <- dat |> subset(dat[[qname]] == 0)
  
  pink <- exp(unqueried[,yname] * (beta0 + beta2 * unqueried[,zname]) - 
                exp(beta0 + beta2 * unqueried[,zname]) -
                (nu0 + nu1 * unqueried[,xstarname] + nu2 * unqueried[,zname])) + 
    exp(unqueried[,yname] * (beta0 + beta1 + beta2 * unqueried[,zname]) -
          exp(beta0 + beta1 + beta2 * unqueried[,zname]))
  
  green <- unqueried[,yname] * (beta0 + beta2 * unqueried[,zname]) - 
    exp(beta0 + beta2 * unqueried[,zname]) -
    (nu0 + nu1 * unqueried[,xstarname] + nu2 * unqueried[,zname])
  
  orange <- unqueried[,yname] * (beta0 + beta1 + beta2 * unqueried[,zname]) -
    exp(beta0 + beta1 + beta2 * unqueried[,zname])
  
  purple <- -1 * (nu0 + nu1 * unqueried[,xstarname] + nu2 * unqueried[,zname])
  
  db0 <- queried[,yname] - exp(beta0 + beta1*queried[,xname] + 
                                 beta2*queried[,zname]) + #queried piece
    (exp(green) * (unqueried[,yname] - exp(beta0 + beta2*unqueried[,zname])) + 
       exp(orange) * (unqueried[,yname] - exp(beta0 + beta1 + beta2*unqueried[,zname])))/
    (pink) |> #unqueried
    sum()
  
  db1 <- queried[,yname] * queried[,xname] - queried[,xname] * 
    exp(beta0 + beta1*queried[,xname] + beta2*queried[,zname]) + #queried piece
    exp(orange) * (unqueried[,yname] - exp(beta0 + beta1 + beta2 * unqueried[,zname]))/ 
    (pink) |> #unqueried
    sum()  
  
  db2 <- queried[,yname] * queried[,zname] - queried[,zname] *
    exp(beta0 + beta1*queried[,xname] + beta2*queried[,zname]) + #queried piece
    exp(green) * (unqueried[,yname] * unqueried[,zname] - 
                      unqueried[,zname] * exp(beta0 + beta2 * unqueried[,zname])) +
    exp(orange) * (unqueried[,yname] * unqueried[,zname] - unqueried[,zname] * 
                       exp(beta0 + beta1 + beta2 * unqueried[,zname])) /
    (pink) |> #unqueried piece
    sum()
  
  dn0 <- (-1 * exp(-1 * (nu0 + nu1*queried[,xstarname] + nu2*queried[,zname]))) / 
    (1 + exp(-1 * (nu0 + nu1*queried[,xstarname] + nu2*queried[,zname]))) -
    (1 - queried[,xname]) + #queried piece
    -1 * exp(green) / 
    pink +
    exp(purple) /
    (1 + exp(purple)) |> #unqueried piece
    sum()
  
  dn1 <- (-1 * queried[,xstarname] * exp(-1 * (nu0 + nu1*queried[,xstarname] + nu2*queried[,zname]))) /
    (1 + exp(-1 * (nu0 + nu1*queried[,xstarname] + nu2*queried[,zname]))) -
    (1 - queried[,xname]) * queried[,xstarname] + #queried piece
    (-1 * unqueried[,xstarname] * exp(green)) /
    pink +
    (unqueried[,xstarname] * exp(purple)) /
    (1 + exp(purple)) |> #unqueried piece
    sum()
  
  dn2 <- (-1 * queried[,zname] * exp(-1 * (nu0 + nu1*queried[,xstarname] + nu2*queried[,zname]))) /
    (1 + exp(-1 * (nu0 + nu1*queried[,xstarname] + nu2*queried[,zname]))) -
    (1 - queried[,xname]) * queried[,zname] + #queried piece
    (-1 * unqueried[,zname] * exp(green)) /
    pink +
    (unqueried[,zname] * exp(purple)) /
    (1 + exp(purple)) |> #unqueried piece
    sum()
  
  return(c(db0, db1, db2, dn0, dn1, dn2))
}