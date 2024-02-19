### PETRELLA - ASSIGNMENT 2 - exercise 3

## (a)
# Set the seed for reproducibility
set.seed(123)

# Define the parameters
sigma_e <- 1
sigma_eta <- 1
sigma_xi <- 1
m <- 1
b <- 1

# Set the length of the time series
n <- 300

# Generate the errors
epsilon <- rnorm(n, 0, sigma_e)
eta <- rnorm(n, 0, sigma_eta)
xi <- rnorm(n, 0, sigma_xi)

ts.plot(epsilon)
lines(eta, col = "red")
lines(xi, col = "blue")

# Generate the time series
mu_t <- 0
beta_t <- 0
y <- rep(NA, n)

for (i in 1:n) {
  y[i] <- mu_t + epsilon[i]
  mu_t <- beta_t + mu_t + eta[i]
  beta_t <- beta_t + xi[i]
}

# Plot the time series
plot(y, type = "l", xlab = "Time", ylab = "y", main = "LLT Model Simulation")


acf(y,lag.max = 100, drop.lag.0 = FALSE)
pacf(y,lag.max = 10, drop.lag.0 = FALSE)

# lag operator
dy_sq = diff(y)^2
ts.plot(dy_sq)
acf(dy_sq,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy_sq,lag.max = 100, drop.lag.0 = FALSE)

##  (b)
# (sigma_e, sigma_eta) = (0.1,1);(1,1);(1,0.1)
# Define the parameters
phi <- 0.8
mu <- 0
sigma_e <- 1
sigma_eta <- 1

# Set the length of the time series
n <- 300

# Generate the errors
epsilon <- rnorm(n, 0, sigma_e)
ts.plot(epsilon)
acf(epsilon, lag.max = 80, drop.lag.0 = FALSE)

eta <- rnorm(n, 0, sigma_eta)
ts.plot(eta)
acf(eta, drop.lag.0 = FALSE, lag.max = 120)

ts.plot(eta)
lines(epsilon, col = "red")

# Initialize the time series
y <- rep(NA, n)

# Set the initial value of mu
mu_t <- rnorm(1, mu, sqrt(sigma_eta / (1 - phi^2)))

# Generate the time series
for (i in 1:n) {
  y[i] <- mu_t + epsilon[i]
  mu_t <- phi * mu_t + eta[i]
}

ts.plot(y)
lines(mu_t,col="red")

acf(y,lag.max = 100, drop.lag.0 = FALSE)
pacf(y,lag.max = 100, drop.lag.0 = FALSE)

