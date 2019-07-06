source("c:\\Users\\wrigh\\OneDrive\\Desktop\\CWB\\CWB2019code.R")

## Question 1

### 1.1

# negloglike: Computes negative log-likelihood
# Input:
# param : vector with 2 values, N and theta
# Y : vector with 2 values, y1 and y2
# Output :
# a scalar equal to negative log-likelihood

negloglike <- function(param, Y) {
  N <- param[1]
  theta <- param[2]
  if (N >= max(Y)) {
    nll <- sum(lgamma((Y + 1))) +
      sum(lgamma(N - Y + 1)) -
      2 * lgamma(N + 1) +
      2 * N * log(1 + exp(theta)) -
      (sum(Y)) * theta
    return(nll)
  } else {
    return(+Inf)
  }
}

### 1.2
# Calculating maximum likelihood estimates for N and theta
Y <- c(256, 237)
opt <- optim(
  par = c(2 * max(Y), logit(0.5)),
  fn = negloglike, Y = Y
)

theta_hat <- opt$`par`[2]
phi <- ilogit(theta_hat)
phi

### 1.3
# Computing Hessian matrix
estimates <- opt$`par`
hessian <- optimHess(
  par = estimates,
  fn = negloglike, Y = Y
)
hessian

# Constructing 95% confidence interval for N
z_0975 <- qnorm(0.975)
inversehessian <- solve(hessian)
se <- sqrt(diag(inversehessian))
se_N <- se[1]
N_hat <- opt$`par`[1]
N <- max(Y) # We estimate this by the prelude as we know there must be at least 256
# at worst this is the upper and lower bound as we are dividing by n

lower_bound <- max(Y)
interval <- c(N_hat - z_0975 * se_N, N_hat + z_0975 * se_N)

# This does respect the lower bound

## Question 2
### 2.1
# See paper and put in markdown

### 2.2
# negloglike: Compute
# Input:
# param : vector with 2 values, N and theta
# Y : vector with 2 values, y1 and y2
# Output :
# a hessian matrix
myhessian <- function(param, Y) {
  N <- param[1]
  theta <- param[2]
  second_partial_theta <- 2 * N * exp(theta) / (1 + exp(theta))^(2)
  second_partial_N <- sum(psigamma(N - Y + 1, 1)) - 2 * psigamma(N + 1, 1)
  partial_theta_N <- 2 * exp(theta) / (1 + exp(theta))
  hess <- matrix(c(
    second_partial_N, partial_theta_N,
    partial_theta_N, second_partial_theta
  ),
  nrow = 2,
  ncol = 2
  )
  return(hess)
}
# Calculating relative difference between elements
myhess <- myhessian(estimates, Y)
relative_difference <- (myhess - hessian) / hessian

### 2.3

abs_diff_L2N <- abs(myhess[1, 1] - hessian[1, 1])
abs_diff_L2N

# DO WE NEED TO CHANGE TO H=0.0001 SOMEWHERE
# Comparing absolute difference in hessian evaluations using bound for second order difference approximations
L0 <- negloglike(estimates, Y)
L1 <- sum(digamma(N_hat - Y + 1)) - 2 * digamma(N_hat + 1) + 2 * log(1 + exp(theta_hat))
L4 <- sum(psigamma(N_hat - Y + 1, 3)) - 2 * psigamma(N_hat + 1, 3)
h <- 0.0001

myhessian(estimates, Y)[1, 1] - hessian[1, 1]

epsilon_0 <- .Machine$double.eps
bound <- epsilon_0 * (4 * L0 + 2 * abs(theta_hat) * L1) / (h^2) + h^2 * L4 / 12
bound

### 2.4

# AVERAGE THIS OVER A NUMBER OF RUNS
mb <- microbenchmark::microbenchmark(myhessian(estimates, Y), optimHess(
  par = estimates,
  fn = negloglike, Y = Y
))
mb <- summary(mb)
mb

## Question 3
### 3.1
opti <- function(Y) {
  MLE <- optim(
    par = c(2 * max(Y), 0),
    fn = negloglike,
    Y = Y
  )$par
  return(MLE)
}

arch_boot <- function(param, J) {
  N <- floor(param[1])
  n <- length(param)
  phi_hat <- ilogit(param[2])
  Y <- rbinom(2 * J, N, phi_hat) #2 pairs of left and right pairs of Y
  samples <- matrix(Y, J, n) # Storing random values in matrix
  matrix <- t(apply(samples, 1, opti)) # Apply optim to every row and transpose matrix
  return(matrix)
}

### 3.2

# Using bootstrap principle
J <- 10000
BS_mles <- arch_boot(estimates, J)
mean_N_hat <- mean(BS_mles[, 1])
bias_N <- mean_N_hat - estimates[1]
sd_N_hat <- sd(BS_mles[, 1] - estimates[1])

mean_theta_hat <- mean(BS_mles[, 2])
bias_theta <- mean_theta_hat - estimates[2]
sd_theta_hat <- sd(BS_mles[, 2] - estimates[2])

### 3.3

# Confidence interval for log N using bootstrap
logdiff_BS_mles <- log(BS_mles[, 1]) - log(estimates[1])
quantile_N <- quantile(logdiff_BS_mles, probs = c(0.975, 0.025))

CI_log_N <- log(estimates[1]) - quantile_N
CI_N <- exp(CI_log_N)

diff_theta <- BS_mles[, 2] - estimates[2]
quantile_theta <- quantile(diff_theta, probs = c(0.975, 0.025))
CI_theta <- estimates[2] - quantile_theta
CI_phi <- ilogit(CI_theta)

#Change interval names to opposite order

## Question 4
# 4.1

# see paper do not LOSE

# 4.2

data <- read_TMIN_data()

J <- 20
SE <- c(1:J)
DS <- c(1:J)
Brier <- c(1:J)

S_boot_train <- data.frame(
  SE = numeric(J),
  DS = numeric(J),
  Brier = numeric(J)
)
S_boot_test <- data.frame(
  SE = numeric(J),
  DS = numeric(J),
  Brier = numeric(J)
)

K <- 10

for (j in 1:J) {
  data_resample <- data_list_resample(data)
  cwb4 <- cwb4_scores(data_resample, K)
  S_boot_train[j, ] <- cwb4$cvk_mean
  S_boot_test[j, ] <- cwb4$test
}

head(S_boot_train)
head(S_boot_test)

# 4.3

train_sample_means_scores <- colMeans(S_boot_train)
train_sample_sd <- apply(S_boot_train, 2, sd)

test_sample_means_scores <- colMeans(S_boot_test)
test_sample_sd <- apply(S_boot_test, 2, sd)

# 4.4

bias_score <- colMeans(S_boot_train - S_boot_test)
standard_devation_score <- apply(S_boot_train-S_boot_test,2,sd)

#Variablity dominates the score estimates as it is one order higher
#Yes the cross validation scores exceed the test scores
