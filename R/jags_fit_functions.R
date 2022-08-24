FitOneVariableParamSpikeAndSlab <- function(data, nu_0) {
  ni <- length(unique(data$i))
  nu_1 <- 1000

  model_code <- "{
    for (i in 1:ni){
      delta[i] ~ dbern(theta)
    }
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })
  "
}
