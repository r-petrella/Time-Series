#### Petrella Assignment 3
# Exercise (1)

################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
################################################################################

# write detailed comments each line of the following code and improve it whenever you
# you feel able to do that 

n = 300    # length of the series

sigma_eta = c(0.1,1)    # variance of the error term which 
# captures the variation in the time series that is not explained by the trend 
# component of the model

eta = sigma_eta[1] * rnorm(n) # to get the normal distribution of the error for the length
# of the time series


sigma_e = c(1, 0.1) # represents the variance of the residual associated with 
# the other explanatory variables in the model



q = sigma_eta[1]^2/sigma_e[1]^2   # signal-to-noise ratio  

e = sigma_e[1] * rnorm(n)   # to generate the norm distribution of the error for the length
# of the time series

#  nu = c(3,  6, 12, 28, 200)
#  e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu[1])
#  ts.plot(e)


ts.plot(e)   # to plot the time series
lines(eta, col = "red")   # to write on the plot above the distribution of eta

phi <- 0.8   # weigth of the coefficient for mu[t]
mu = 0  # the initial value of mu
mu[1] <- 0    # to create a list for mu
y = 0 # initial value of y
for(t in 1:(n-1)){ # we create a for loop of t from 1 to 299
  mu[t+1] = phi * mu[t] + eta[t] # AR(1) process (transition equation)
  y[t+1] =  mu[t+1] + e[t+1]   # observed equation
}
ts.plot(mu)   # plot of mu

# y = mu + e

ts.plot(y)  # plot of the observations of the process
lines(mu,col="red") # representation of mu over the above plot


acf(y, lag.max = 60, drop.lag.0 = FALSE)  # auto-correlation function of y (max 60 points)
pacf(y, lag.max = 60)  # partial auto-correlation function (max 60 points)
acf(mu, lag.max = 60, drop.lag.0 = FALSE) # auto-correlation function of mu
pacf(mu, lag.max = 60)   # partial auto-correlation function


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################


mu_pred <- 0   # initial predicted state estimate. It is an estimate of the 
# underlying state of the system given the observed data up to time t

P <- 0      # The state error covariance matrix. This matrix represents the 
# uncertainty in the current state estimate. At each time step, the Kalman 
# filter updates P based on the predicted and observed data

v <-0     # The innovation or measurement residual. It represents the  
#difference between the observed measurement and the predicted measurement
# based on the current state estimate. The Kalman filter uses v 
# to update the state estimate

K <- 0    # The Kalman gain. It is a weighting factor that determines how
# much to trust the predicted state estimate versus the observed measurement.
# The Kalman filter uses K to update the state estimate based on v     

F <- 0  # F: The innovation covariance matrix. It represents the uncertainty in
# the observed measurement given the current state estimate. The Kalman filter
# uses F to calculate K and update the state estimate.           

llk<-0  # The log-likelihood of the observed data given the state estimate.
# The Kalman filter uses llk to assess the goodness of fit of the filter.
# The log-likelihood is updated at each time step based on the predicted
# and observed data.     


mu_pred[1] = 0   # use mu_pred as list
P[1] = (sigma_eta[1]^2)/(1-phi^2)  # variance of AR(1) process ([1,1] element of P)
llk[1]=0  # initial log-likelihood


for(t in 1:(n-1)){   # creating a for loop to implement the kalman filter for a 
  # linear state-space model
  v[t] = y[t] - mu_pred[t]; # calculate the innovation or measurement residual
  F[t] = P[t] + sigma_e[1]^2; #  calculate the innovation covariance matrix
  K[t] = (phi * P[t])/F[t] # calculate the Kalman gain
  P[t+1] = phi^2 * P[t] + sigma_eta[1]^2 - K[t]*F[t]*K[t] # update the error 
  # covariance matrix of the predicted state estimate
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t] # update the predicted state estimate
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t])) # update the log-likelihood of the 
  # observed data
}

llk = sum(llk)  # total likelihood

ts.plot(y)  # plot of the process
lines(mu, col = 'red')  # line of mu over y
lines(mu_pred, col = 'green') #line of the predicted state estimate over y


################################ Exercise (2)
library(MASS)

# Define LGSSM parameters
F <- matrix(c(0.8, 0.2, 0, 1), nrow = 2) # state transition matrix
H <- matrix(c(1, 0), nrow = 1) # observation matrix
Q <- matrix(c(0.1, 0, 0, 0.1), nrow = 2) # state noise covariance matrix
R <- matrix(0.5, nrow = 1) # observation noise covariance matrix
mu_0 <- c(0, 0)  # initial state mean vector
Sigma_0 <- matrix(c(1, 0, 0, 1), nrow = 2) # initial state cov matrix

# Simulate data for LGSSM
set.seed(123)
n_obs <- 100
state <- mvrnorm(n_obs + 1, mu_0, Sigma_0)
obs <- matrix(0, nrow = n_obs, ncol = 1)
for (i in 1:n_obs) {
  obs[i,] <- rnorm(1, H %*% state[i + 1,], sqrt(R))
}

# Implement Kalman filter to estimate predicted state
alpha <- matrix(0, nrow = n_obs + 1, ncol = 2)
V <- list()
alpha[1,] <- mu_0
V[[1]] <- Sigma_0
for (i in 1:n_obs) {
  alpha[i + 1,] <- F %*% alpha[i,]
  V[[i + 1]] <- F %*% V[[i]] %*% t(F) + Q
  K <- V[[i + 1]] %*% t(H) %*% solve(H %*% V[[i + 1]] %*% t(H) + R)
  alpha[i + 1,] <- alpha[i + 1,] + K %*% (obs[i,] - H %*% alpha[i + 1,])
  V[[i + 1]] <- (matrix(1, nrow = 2, ncol = 2) - K %*% H) %*% V[[i + 1]]
}
ts.plot(obs)
# Print estimated predicted state
cat("Predicted state for t+1|t:", alpha[n_obs + 1,], "\n")

