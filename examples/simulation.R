source("corcig.R")
library("vars")

set.seed(234)

# ------------------
# 1. Data Processing
# ------------------

# Number of time series observations
tt <- 1000

# Structural coefficients
Phi <- diag(1, 2)
Phi[lower.tri(Phi)] <- c(0.5)

# Coefficient matrices
Phi1 <- matrix(
  c(0.5, -0.2, 0, 0.6), nrow=2, ncol=2
)
Phi2 <- matrix(
  c(0.25, 0.1, 0.5, 0), nrow=2, ncol=2
)

# Structural matrices
Phi1.star = solve(Phi)%*%Phi1
Phi2.star = solve(Phi)%*%Phi2

# Generate series
series <- matrix(0, 2, tt + 1) 
for (i in 3:(tt + 1)){
  series[, i] <- (
    Phi1.star%*%series[, i-1] + 
      Phi2.star%*%series[, i-2] +  
      solve(Phi)%*%rnorm(2, 0, 1)
  )
}
series <- tail(series, n=tt)
series <- ts(t(series)) # Convert to time series object
dimnames(series)[[2]] <- c("y", "x") # Rename variables

# Plot the series
plot.ts(
  series[(tt-100):tt,], 
  ylab="Process Value", 
  main = "Simulated Time Series", 
  plot.type = "single", 
  col = 1:2
)
legend("topleft",c("y.t","x.t"),cex=.8,col=c("red","black"), pch=c(16,16))

# Estimate VAR(2)
var_est <- VAR(series, p = 2, type = "none")
var_est

# Estimate SVAR(2)
B <- diag(1, 2)
B[lower.tri(B)] <- NA
svar_est <- SVAR(var_est, Bmat = B, max.iter = 1000)
solve(svar_est$B)

# Build design matrix
X.t = cbind(series[,"y"], series[,"x"])
X = X.t
for (i in 1:3) {
  temp = stats::lag(X.t, i)
  X = cbind(temp, X)
}
colnames(X) <- c("y.t", "x.t", "y.t1", "x.t1", "y.t2", "x.t2", "y.t3", "x.t3")
X <- X[complete.cases(X),]

# Check lagged variables
tail(X)

# ------------------
# 2. Calculating CIG
# ------------------

# estimate the covariance matrix
S = cov(X)
# get significance theshold
crit.r = critical.r(nrow(X), alpha=0.05)
# build CIG from covariance matrix
g <- cor2cig(S, crit.r)

# plot full CIG (not just links to variables at time t)
plot(g, layout=layout_with_kk, label.cex=.25, label.dist=10)

# TODO: 
