### Problem 1
df <- MPV::p12.8
n <- 5
set.seed(1); chosen_rows <- sort(sample(seq(1, nrow(df)), n))
df <- df[chosen_rows,]

# part (a)
plot(df$x, df$y)
f_0 <- function(x, theta1, theta2) {
  return(theta1 * exp(theta2 * x))
}
x_seq <- seq(0, 10, length.out = 1e3)
lines(x_seq, f_0(x = x_seq, theta1 = 1, theta2 = 0.5))

# Functions start
theta1_deriv_calc <- function(x, theta2) {
  return(exp(theta2 * x))
}
theta2_deriv_calc <- function(x, theta1, theta2) {
  return(theta1 * x * exp(theta2 * x))
}
cutoff_fun <- function(theta_new, theta_old, delta) {
  cutoff_calc <- (theta_new - theta_old)  / theta_hat_old
  
  # Cutoff when both less than delta
  if (cutoff_calc[1] < delta & cutoff_calc[2] < delta) {
    return(TRUE)
  } else {
    return(c(FALSE, round(cutoff_calc[1], 4), round(cutoff_calc[2], 4)))
  }
}
beta_hat_calc <- function(x, y, Z, theta1, theta2) {
  return(solve(t(Z) %*% Z) %*% t(Z) %*%
           (y - f_0(x = x, theta1 = theta1, theta2 = theta2)))
}
RSS_calc <- function(y, y_hat) {
  return(sum((y_hat - y)^2))
}
# Functions end

# Initialize theta-0
theta1 <- 1; theta2 <- 0.5
theta_vec <- c(theta1, theta2)

# Initialize derivatives, Z
theta1_deriv <- theta1_deriv_calc(x = df$x, theta2 = theta_vec[2])
theta2_deriv <- theta2_deriv_calc(x = df$x,
  theta1 = theta_vec[1], theta2 = theta_vec[2])
Z <- cbind(theta1_deriv, theta2_deriv)

# Initialize beta-hat
beta_hat <- solve(t(Z) %*% Z) %*% t(Z) %*% df$y

# Generate theta-hat-1
theta_hat_update <- beta_hat + theta_vec
theta_hat_old <- theta_vec

# Initialize cutoff
delta <- 1e-6
cutoff <- cutoff_fun(theta_new = theta_hat_update,
                     theta_old = theta_hat_old,
                     delta = delta)

# Initial RSS
y_hat <- f_0(x = df$x,
             theta1 = theta_hat_update[1], theta2 = theta_hat_update[2])
RSS <- RSS_calc(y = df$y, y_hat = y_hat)

counter <- 0
while(!cutoff[1]) {
  # Update counter
  counter <- counter + 1
  print(paste0('Counter: ', counter))
  # Calculate Z
  theta1_deriv <- theta1_deriv_calc(x = df$x, theta2 = theta_hat_update[2])
  theta2_deriv <- theta2_deriv_calc(x = df$x,
    theta1 = theta_hat_update[1], theta2 = theta_hat_update[2])
  Z <- cbind(theta1_deriv, theta2_deriv)
  
  # Calculate beta-hat
  beta_hat <- beta_hat_calc(x = df$x, y = df$y, Z = Z,
    theta1 = theta_hat_update[1], theta2 = theta_hat_update[2])
  
  # theta_{k+1} = theta_k + beta_k
  theta_hat_old <- theta_hat_update
  theta_hat_update <- beta_hat + theta_hat_update
  
  # Calculate cutoff
  cutoff <- cutoff_fun(theta_new = theta_hat_update, theta_old = theta_hat_old, delta = delta)
  print(paste0('Cutoff:', cutoff[2], '; ', cutoff[3]))
  
  # Print update values
  print(paste0('theta1:',
               round(theta_hat_update[1], 4),
               '; theta2',
               round(theta_hat_update[2], 4)))
  
  # Print RSS
  y_hat <- f_0(x = df$x,
    theta1 = theta_hat_update[1], theta2 = theta_hat_update[2])
  RSS <- RSS_calc(y = df$y, y_hat = y_hat)
  print(paste0('RSS:', round(RSS, 4)))
}

plot(df$x, df$y,
     main = 'Initial estimate and updated estimate of parameters',
     xlab = 'x', ylab = 'y')
legend("topleft", legend = c('Initial', 'Update'),
       lty = c(2,2), col = c('black', 'red'))
lines(x_seq, f_0(x = x_seq, theta1 = 1, theta2 = 0.5), lty = 2)
lines(x_seq, f_0(x = x_seq,
  theta1 = theta_hat_update[1],
  theta2 = theta_hat_update[2]),
  col = 'red', lty = 2)

# part (b)
SS_T <- t(df$y) %*% df$y - ((sum(df$y))^2) / n
SS_model <- SS_T - RSS
MS_Res <- RSS / (n - 2)

F_0 <- (SS_model / 2) / MS_Res
qf(0.95, df1 = 2, df2 = 3)

# part (c)
cov_mat <- MS_Res * solve(t(Z) %*% Z)
sqrt(diag(cov_mat))

# part (d)
C <- solve(t(Z) %*% Z)
standard_errors <- sqrt(diag(MS_Res * C))
t_test_statics <- theta_hat_update / standard_errors
alpha <- 0.05
qt(1 - alpha/2, df = 3)

# part (e)
e <- df$y - y_hat
beta_hat_calc <- function(X, y) {
  X <- as.matrix(X)
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
  return(beta_hat)
}
norm_prob_plot <- function(residual_var, x_label,
                           main_title = 'Normal Probability Plot',
                           y_label = 'Probability', n_size=n) {
  ones <- rep(1, n)
  sorted_residuals <- sort(residual_var)
  cumulative_probability <- (1:n_size - 0.5) / n_size
  plot(sorted_residuals, cumulative_probability, main = main_title,
       xlab = x_label,
       ylab = y_label)
  X_temp <- cbind(ones, sorted_residuals)
  beta_hat_temp <- beta_hat_calc(X=X_temp,y=cumulative_probability)
  abline(beta_hat_temp)
}
norm_prob_plot(residual_var = e, x_label = 'Sorted Residuals')
order(e, decreasing = FALSE)
e[order(e, decreasing = FALSE)]

res_vs_fitted_plot <- function(residual_var,
                               main_title,
                               y_label,
                               x_label = 'Predicted Response',
                               pred_response = y_hat) {
  plot(pred_response, residual_var, main = main_title,
       xlab = x_label,
       ylab = y_label,
       ylim = c(min(residual_var)-sd(residual_var),
                max(residual_var)+sd(residual_var)))
}
res_vs_fitted_plot(residual_var = e,
                   main_title = 'Residuals vs. Predicted Response',
                   y_label = 'Residuals')


### Problem 2
df <- MPV::p12.11
n <- 14
set.seed(1); chosen_rows <- sort(sample(seq(1, nrow(df)), n))
df <- df[chosen_rows,]

# part (a)
plot(df$x, df$y, main = 'Problem 12.11 Y vs. X', xlab = 'x', ylab = 'y')

# part (b)
f_0 <- function(x, theta1, theta2, theta3) {
  return(theta1 - theta2 * exp(- theta3 * x))
}

x_seq <- seq(0, 40, length.out = 1e3)
lines(x_seq, f_0(x = x_seq, theta1 = 0.4, theta2 = -0.3, theta3 = 0.15))

# Functions start
theta1_deriv_calc <- function(x) {
  return(rep(1, length(x)))
}
theta2_deriv_calc <- function(x, theta3) {
  return(-exp(-theta3 * x))
}
theta3_deriv_calc <- function(x, theta2, theta3) {
  return(theta2 * x * exp(-theta3 * x))
}
cutoff_fun <- function(theta_new, theta_old, delta) {
  cutoff_calc <- (theta_new - theta_old)  / theta_hat_old
  
  # Cutoff when both less than delta
  if (cutoff_calc[1] < delta &
      cutoff_calc[2] < delta &
      cutoff_calc[3] < delta) {
    return(TRUE)
  } else {
    return(c(FALSE,
             round(cutoff_calc[1], 4),
             round(cutoff_calc[2], 4),
             round(cutoff_calc[3], 4)))
  }
}
beta_hat_calc <- function(x, y, Z, theta1, theta2, theta3) {
  return(solve(t(Z) %*% Z) %*% t(Z) %*%
           (y - f_0(x = x, theta1 = theta1, theta2 = theta2, theta3 = theta3)))
}
RSS_calc <- function(y, y_hat) {
  return(sum((y_hat - y)^2))
}
# Functions end

# Initialize theta-0
theta1 <- 0.5; theta2 <- -0.3; theta3 <- 0.15
theta_vec <- c(theta1, theta2, theta3)

# Initialize derivatives, Z
theta1_deriv <- theta1_deriv_calc(x = df$x)
theta2_deriv <- theta2_deriv_calc(x = df$x, theta3 = theta_vec[3])
theta3_deriv <- theta3_deriv_calc(x = df$x,
                                  theta2 = theta_vec[2], theta3 = theta_vec[3])

Z <- cbind(theta1_deriv, theta2_deriv, theta3_deriv)

# Initialize beta-hat
beta_hat <- solve(t(Z) %*% Z) %*% t(Z) %*% df$y

# Generate theta-hat-1
theta_hat_update <- beta_hat + theta_vec
theta_hat_old <- theta_vec

# Initialize cutoff
delta <- 1e-6
cutoff <- cutoff_fun(theta_new = theta_hat_update,
                     theta_old = theta_hat_old,
                     delta = delta)

# Initial RSS
y_hat <- f_0(x = df$x,
             theta1 = theta_hat_update[1],
             theta2 = theta_hat_update[2],
             theta3 = theta_hat_update[3])
RSS <- RSS_calc(y = df$y, y_hat = y_hat)

# Run the while-loop
counter <- 0
while(!cutoff[1]) {
  # Update counter
  counter <- counter + 1
  print(paste0('Counter: ', counter))
  # Calculate Z
  theta1_deriv <- theta1_deriv_calc(x = df$x)
  theta2_deriv <- theta2_deriv_calc(x = df$x, theta3 = theta_hat_update[3])
  theta3_deriv <- theta3_deriv_calc(x = df$x,
                                    theta2 = theta_hat_update[2], theta3 = theta_hat_update[3])
  Z <- cbind(theta1_deriv, theta2_deriv, theta3_deriv)
  
  # Calculate beta-hat
  beta_hat <- beta_hat_calc(x = df$x, y = df$y, Z = Z,
                            theta1 = theta_hat_update[1],
                            theta2 = theta_hat_update[2],
                            theta3 = theta_hat_update[3])
  
  # theta_{k+1} = theta_k + beta_k
  theta_hat_old <- theta_hat_update
  theta_hat_update <- beta_hat + theta_hat_update
  
  # Calculate cutoff
  cutoff <- cutoff_fun(theta_new = theta_hat_update, theta_old = theta_hat_old, delta = delta)
  print(paste0('Cutoff:', cutoff[2], '; ', cutoff[3], '; ', cutoff[3]))
  
  # Print update values
  print(paste0('theta1:', round(theta_hat_update[1], 4),
               '; theta2', round(theta_hat_update[2], 4),
               '; theta2', round(theta_hat_update[2], 4)))
  
  # Print RSS
  y_hat <- f_0(x = df$x,
               theta1 = theta_hat_update[1],
               theta2 = theta_hat_update[2],
               theta3 = theta_hat_update[3])
  RSS <- RSS_calc(y = df$y, y_hat = y_hat)
  print(paste0('RSS:', round(RSS, 4)))
}

plot(df$x, df$y,
     main = 'Initial estimate and updated estimate of parameters',
     xlab = 'x', ylab = 'y')
legend("topright", legend = c('Initial', 'Update'),
       lty = c(2,2), col = c('black', 'red'))
x_seq <- seq(0, 40, length.out = 1e3)
lines(x_seq, f_0(x = x_seq, theta1 = 0.4, theta2 = -0.3, theta3 = 0.15), lty = 2)
lines(x_seq, f_0(x = x_seq,
                 theta1 = theta_hat_update[1],
                 theta2 = theta_hat_update[2],
                 theta3 = theta_hat_update[3]), lty = 2, col = 'red')

# part (c)
SS_T <- t(df$y) %*% df$y - ((sum(df$y))^2) / n
SS_model <- SS_T - RSS
MS_Res <- RSS / (n - 3)

F_0 <- (SS_model / 3) / MS_Res
qf(0.95, df1 = 3, df2 = n - 3)

# part (d)
C <- solve(t(Z) %*% Z)
standard_errors <- sqrt(diag(MS_Res * C))
alpha <- 0.05
t_value <- qt(1 - alpha / 2, n - 3)
round(theta_hat_update + t_value * standard_errors, 4)
round(theta_hat_update - t_value * standard_errors, 4)

# part (e)
e <- df$y - y_hat
beta_hat_calc <- function(X, y) {
  X <- as.matrix(X)
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
  return(beta_hat)
}
norm_prob_plot <- function(residual_var, x_label,
                           main_title = 'Normal Probability Plot',
                           y_label = 'Probability', n_size=n) {
  ones <- rep(1, n)
  sorted_residuals <- sort(residual_var)
  cumulative_probability <- (1:n_size - 0.5) / n_size
  plot(sorted_residuals, cumulative_probability, main = main_title,
       xlab = x_label,
       ylab = y_label)
  X_temp <- cbind(ones, sorted_residuals)
  beta_hat_temp <- beta_hat_calc(X=X_temp,y=cumulative_probability)
  abline(beta_hat_temp)
}
norm_prob_plot(residual_var = e, x_label = 'Sorted Residuals')
order(e, decreasing = FALSE)
e[order(e, decreasing = FALSE)]

res_vs_fitted_plot <- function(residual_var,
                               main_title,
                               y_label,
                               x_label = 'Predicted Response',
                               pred_response = y_hat) {
  plot(pred_response, residual_var, main = main_title,
       xlab = x_label,
       ylab = y_label,
       ylim = c(min(residual_var)-sd(residual_var),
                max(residual_var)+sd(residual_var)))
}
res_vs_fitted_plot(residual_var = e,
                   main_title = 'Residuals vs. Predicted Response',
                   y_label = 'Residuals')
