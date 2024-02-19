###############
## ASSIGNMENT N PETRELLA RICCARDO
###############


## EXERCISE 2

### KALMAN FILTER WITH PARAMETER OPTIMIZATION ###
library("TSA")
data("Nile")
ts.plot(Nile)
y = Nile
ts.plot(y)
y_index = time(y)

# initiating variables for Kalman filter algorithm
mu_pred <- 0   
P <- 0         
v <-0          
K <- 0         
F <- 0           
llk<-0 

sigma_e = sigma_eta = sd(y)

LLM_kalman_filter = function(y, sigma_EPS = sigma_e, sigma_ETA = sigma_eta){
  # setting initial conditions for predicted hidden state, predicted variance of 
  # hidden state, and log-likelihood
  n = length(y)
  mu_pred[1] = mean(y)
  P[1] = var(y)
  llk <- 0
  dllk <- numeric(length = n)
  
  # implementing Kalman filter recursion
  for(t in 1:(n-1)){
    # computing prediction error
    v[t] = y[t] - mu_pred[t]; 
    # computing predicted variance of prediction error
    F[t] = P[t] + sigma_EPS^2;
    # computing Kalman gain
    K[t] = P[t]/F[t]
    # updating predicted variance of hidden state
    P[t+1] = P[t] + sigma_ETA^2 - F[t]*(K[t])^2
    # updating predicted hidden state
    mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
    # computing log-likelihood
    dllk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
    llk  = llk + dllk[t] ### ok
  }
  
  return(list(mu_pred = mu_pred, llk = llk, dllk = dllk, v = v, 
              P = P, F = F))
  
}

# defining the negative log-likelihood function to be minimized
neg_log_likelihood <- function(params) {
  s_eps <- params[1]
  s_eta <- params[2]
  
  result <- LLM_kalman_filter(y, s_eps, s_eta)
  return(-result$llk)  # Negative log-likelihood
}

# setting initial parameter values for optimization
initial_params <- c(180, 70)

# performing optimization using the "optim" function
# i tried different combinations for the starting points. These are the ones
# that give the closest result to the ones in figure 2.1 and also the  
# likelihood is close enough to the one in table 2.1
optim_result <- optim(
  par = initial_params,
  fn = neg_log_likelihood,
  method = "BFGS"
)

# extracting the optimized parameter values
optim_par <- optim_result$par
optim_par


# fitting the Kalman filter
kf <- LLM_kalman_filter(y, optim_par[1], optim_par[2])
kf$llk
y_pred <- kf$mu_pred
y_pred = ts(y_pred, start = y_index[1])


# plotting
ts.plot(y, lwd = 2, type = "p")
lines(y_pred, lwd = 2, col = "red")

# plotting the filtered state variance
state_var = ts(kf$P, start = y_index[1])
ts.plot(state_var)

# plotting prediction errors
pred_err = ts(kf$v, start = y_index[1])
ts.plot(pred_err)

# plotting prediction variance
pred_var = ts(kf$F, start = y_index[1])
ts.plot(pred_var)

# the plots seem to be behave as the graphs shown in the textbook. The filtered
# state variance and the prediction variance both stabilize after a few 
# iterations. The residuals too seem to exhibit the same pattern shown by the 
# graphs in the textbook

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## EXERCISE N 3
library("TSA")
data("Nile")
y = Nile
y_index = time(y)

LLM_kalman_filter_smoother <- function(y, sigma_EPS, sigma_ETA){
  n <- length(y)
  
  # Initialize
  mu_pred <- mu_smooth <- P <- P_smooth <- v <- K <- F <- L <- numeric(n)
  mu_pred[1] <- mean(y)
  P[1] <- var(y)
  dllk = rep(0, times = n)
  llk = 0
  
  # Kalman filter
  for(t in 1:(n - 1)) {
    v[t] <- y[t] - mu_pred[t] # prediction error
    F[t] <- P[t] + sigma_EPS^2 # variance of prediction error
    K[t] <- P[t] / F[t] # Kalman gain
    mu_pred[t + 1] <- mu_pred[t] + K[t] * v[t] # update state prediction
    P[t + 1] <- (1 - K[t]) * P[t] + sigma_ETA^2 # update variance prediction
    L[t] <- 1 - K[t] # learning coefficient
    
    # computing log-likelihood
    dllk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
    llk  = llk + dllk[t] ### ok
  }
  
  # Final prediction error and variance for the filter stage
  v[n] <- y[n] - mu_pred[n]
  F[n] <- P[n] + sigma_EPS^2
  
  # Kalman smoother
  r <- N <- numeric(n)
  mu_smooth[n] <- mu_pred[n]
  P_smooth[n] <- P[n]
  
  for(t in n:2){
    r[t - 1] <- v[t]/F[t] + L[t]*r[t]
    N[t - 1] <- 1/F[t] + L[t]^2 * N[t]
    P_smooth[t] <- P[t] - P[t]^2 * N[t - 1]
    mu_smooth[t] <- mu_pred[t] + P[t] * r[t - 1]
  }
  
  P_smooth[1] <- P[1]
  mu_smooth[1] <- mu_pred[1]
  
  return(list(mu_pred = mu_pred, mu_smooth = mu_smooth, P = P, P_smooth = P_smooth,
              r = r, N = N, llk = llk, F = F, K = K, v = v))
}


# defining the negative log-likelihood function to be minimized
neg_log_likelihood <- function(params) {
  s_eps <- params[1]
  s_eta <- params[2]
  
  result <- LLM_kalman_filter_smoother(y, s_eps, s_eta)
  return(-result$llk)  # Negative log-likelihood
}

# setting initial parameter values for optimization using observed standard 
# deviation. I tried several starting points, 180 and 70 are the ones that 
# yield the closest results in terms of filtering and smooth to the ones in 
# the textbook
initial_params <- c(180, 70)

# performing optimization using the "optim" function
optim_result <- optim(
  par = initial_params,
  fn = neg_log_likelihood,
  method = "BFGS"
)

# extracting the optimized parameter values
optim_par <- optim_result$par
optim_par

kf_smooth = LLM_kalman_filter_smoother(y, abs(optim_par[1]), abs(optim_par[2]))
kf_smooth$llk

m_sm = kf_smooth$mu_smooth
m_sm = ts(m_sm, start = y_index[1])
m_pr = kf_smooth$mu_pred
m_pr = ts(m_pr, start = y_index[1])

P_sm = kf_smooth$P_smooth
P_sm = ts(P_sm, start = y_index[1])
r_sm = kf_smooth$r
r_sm = ts(r_sm, start = y_index[1])
N_sm = kf_smooth$N
N_sm = ts(N_sm, start = y_index[1])
F_ = kf_smooth$F
F_ = ts(F_, start = y_index[1])
v = kf_smooth$v
v = ts(v, start = y_index[1])
K = kf_smooth$K
K = ts(K, start = y_index[1])


## plotting filtered and smoothed state
ts.plot(y, type = "l")
lines(m_pr, lwd = 2, col = "red")
lines(m_sm, lwd = 2, col = "blue")
legend("topright", legend = c("filtered", "smoothed"), col = c("red", "blue"), lwd = 2)

# the filter and the smoother seem to behave the same way as in the textbook's 
# figure 2.2




## SECTION 2.6 SIMULATION
n = length(y)
set.seed(123)
# drawing random deviates for epsilon and eta
e_plus = rnorm(n, 0, abs(optim_par[1]))
eta_plus = rnorm(n, 0, abs(optim_par[2]))

# generating an LLM process using the random deviates
mu = numeric(n)
llm = numeric(n)

# initializing hidden state and observed data
mu[1] <- m_sm[2]
llm[1] <- mu[1] + e_plus[1]

for(t in 1:(n - 1)){
  # updating hidden state
  mu[t+1] = mu[t] + eta_plus[t]
  # updating observed data
  llm[t+1] =  mu[t+1] + e_plus[t+1]
}

# turning into time series object
llm = ts(llm, start = y_index[1])
ts.plot(llm)


# running kalman smoother on simulated data
kf_smooth_sim = LLM_kalman_filter_smoother(llm, abs(optim_par[1]), abs(optim_par[2]))
v_sim = kf_smooth_sim$v
F_sim = kf_smooth_sim$F
r_sm_sim = kf_smooth_sim$r
K_sim = kf_smooth_sim$K



# computing the smoothing error
u_nile = rep(0, times = n)
for (t in 1:n) {
  u_nile[t] = v[t]/F_[t] - K[t]*r_sm[t]
}

u_sim = rep(0, times = n)
for (t in 1:n) {
  u_sim[t] = v_sim[t]/F_sim[t] - K_sim[t]*r_sm_sim[t]
}

# computing state error
e_hat = (optim_par[1]^2)*u_nile
e_hat = ts(e_hat, start = y_index[1])
e_hat_plus = (optim_par[1]^2)*u_sim
e_hat_plut = ts(e_hat_plus, start = y_index[1])

eta_hat = rep(0, times = n)
eta_hat = ts(eta_hat, start = y_index[1])
for (t in 1:n) {
  eta_hat[t] = m_sm[t+1] - m_sm[t]
}  


# computing conditional draw for epsilon_t given Y_t
e_tilde = e_plus - e_hat_plus + e_hat

# computing conditional state given Y_t
alpha_tilde = y - e_tilde

# computing eta_tilde
eta_tilde = rep(0, times = n)
eta_tilde = ts(eta_tilde, start = y_index[1])
for (t in 1:n) {
  eta_tilde[t] = alpha_tilde[t+1] - alpha_tilde[t]
}


### PLOTTING ##
## smoothed state alpha_hat and simulated sample state
# turning mu into a ts object
mu = ts(mu, start = y_index[1])
msm_range <- range(m_sm)  
msm_range[2] <- msm_range[2] * 1.05  
ts.plot(mu, lwd = 2, type = "l", ylim = msm_range)
lines(m_sm, lwd = 2, col = "red")

## smoothed state and sample alpha_tilde
ts.plot(alpha_tilde, type = "l")
lines(m_sm, lwd = 2, col = "red")

## smoothed observation error e_hat and sample e_tilde
ts.plot(e_tilde, type = "l")
lines(e_hat, lwd = 2, col = "red")

## smoothed state error eta_hat and sample eta_tilde
ts.plot(eta_tilde, type = "l")
lines(eta_hat, lwd = 2, col = "red")


# except for the fact that i don't know with which seed the LLM process was
# generated, the results from the plots seem to get along with the graphs in
# figure 2.4 from the textbook

