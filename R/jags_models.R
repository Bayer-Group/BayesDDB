# JAGS Models

# Bayesian model for the calculation of predictive probability of observing
# current mean given external data
model_string_PPP <- "
  model
  {
    for (i in 1:n){
      y[i] ~ dnorm(m,pre)
   }
  m ~ dnorm(0,tau)
  pre ~ dgamma(alpha,beta)
  new ~ dnorm(m,pre*N1)
  }
  "

# Normal model without any discounting, used for the current-trial treatment arm.
model_string_TREATMENT <- "
  model
  {
    for (i in 1:n){
      y[i] ~ dnorm(mu,pre)
    }
    mu ~ dnorm(0,tau)
    pre ~ dgamma(alpha,beta)
  }
  "


# Power-prior normal model implemented with the zeros trick, so that each
# observation's log-likelihood can be weighted by gamma[i].
model_string_CONTROL <- "
  model
  {
    for(i in 1:n){
      # Define minus log-likelihood plus a large constant to ensure it remains positive
      log_lik[i] <- (log(pre)/2 -log(6.283185)/2 -(y[i]-mu)^2*pre/2)*gamma[i]
      phi[i] <- -log_lik[i] + 1000000 
      zeros[i] ~ dpois(phi[i])
    }

  # Prior distributions
    mu ~ dnorm(0, tau)
    pre ~ dgamma(alpha, beta)
  }
  "
