library(MPV)

### Problem 1
n <- 20
set.seed(1); chosen_rows <- sort(sample(seq(1,27), n))
df <- MPV::table.b5
df <- df[chosen_rows,c(1, 2, 7)]
ones <- rep(1, n)
y <- df[,1]
X <- cbind(ones, df[,c(2,3)])

beta_hat_calc <- function(X, y) {
  X <- as.matrix(X)
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
  return(beta_hat)
}
H_calc <- function(X) {
  X <- as.matrix(X)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  return(H)
}
y_hat_calc <- function(H, y) {
  y_hat <- H %*% y
  return(y_hat)
}
e_calc <- function(y, y_hat) {
  e <- y - y_hat
  return(e)
}

beta_hat <- beta_hat_calc(X=X, y=y)
H <- H_calc(X)
y_hat <- y_hat_calc(H=H, y=y)
e <- e_calc(y=y, y_hat=y_hat)

#### part (a)
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

### part (b)
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

### part (c)
# studentized residuals
SS_Res_calc <- function(y, beta_hat, X) {
  SS_Res <- (t(y) %*% y) - (t(beta_hat) %*% t(X) %*% y)
  return(SS_Res)
}
SS_Res <- SS_Res_calc(y = y, beta_hat = beta_hat, X = X)
p <- ncol(X)
MS_Res <- SS_Res / (n - p)
H_diag <- diag(H)
studentized_residuals <- sapply(1:n, function(x) {
  e[x] / sqrt(MS_Res * (1 - H_diag[x]))
})

# res vs. fitted, norm prob res, res vs. x1, res vs. x6
res_vs_regressor <- function(residual_var,
                             main_title1, main_title2,
                             ylabel,
                             X_df=X) {
  x1 <- X_df[,2]; x6 <- X_df[,3]
  plot(x1, residual_var,
       main = main_title1,
       ylab = ylabel,
       xlab = 'Space time, min.')
  abline(h = 0)

  plot(x6, residual_var,
       main = main_title2,
       ylab = ylabel,
       xlab = 'Solvent total')
  abline(h = 0)
}

par(mfrow = c(2,2))
res_vs_fitted_plot(residual_var = studentized_residuals,
                   main_title = 'Studentized Residuals vs. Predicted Response',
                   y_label = 'Studentized Residuals')
norm_prob_plot(residual_var = studentized_residuals, x_label = 'Studentized Residuals')
res_vs_regressor(residual_var = studentized_residuals,
                 main_title1 = 'Studentized Residuals vs. Space time, min',
                 main_title2 = 'Studentized Residuals vs. Solvent total',
                 ylabel = 'Studentized Residuals')

# Analyze the outliers
head(order(x6, decreasing = TRUE), 5)
studentized_residuals[12] # bottom right
studentized_residuals[1]
studentized_residuals[18]
studentized_residuals[8] # top right
studentized_residuals[11]

# R-student
S_squared <- sapply(1:n, function(x) {
  ((n - p) * MS_Res - ((e[x]^2) / (1 - H_diag[x]))) / (n - p - 1)
})

R_student_res <- sapply(1:n, function(x) {
  e[x] / sqrt(S_squared[x] * (1 - H_diag[x]))
})

par(mfrow = c(2,2))
res_vs_fitted_plot(residual_var = R_student_res,
                   main_title = 'R-student vs. Predicted Response',
                   y_label = 'R-student')
norm_prob_plot(residual_var = R_student_res, x_label = 'R-student')
res_vs_regressor(residual_var = R_student_res,
                 main_title1 = 'R-student vs. Space time, min',
                 main_title2 = 'R-student vs. Solvent total',
                 ylabel = 'R-student')

### part (d)
# standardized residuals
standardized_res <- sapply(1:n, function(x) e[x] / sqrt(MS_Res))

par(mfrow = c(2,2))
res_vs_fitted_plot(residual_var = standardized_res,
                   main_title = 'Standardized Residuals vs. Predicted Response',
                   y_label = 'Standardized Residuals')
norm_prob_plot(residual_var = standardized_res, x_label = 'Standardized Residuals')
res_vs_regressor(residual_var = standardized_res,
                 main_title1 = 'Standardized Residuals vs. Space time, min',
                 main_title2 = 'Standardized Residuals vs. Solvent total',
                 ylabel = 'Standardized Residuals')


# PRESS residuals
PRESS_res <- sapply(1:n, function(x) e[x] / (1 - H_diag[x]))

par(mfrow = c(2,2))
res_vs_fitted_plot(residual_var = PRESS_res,
                   main_title = 'PRESS Residuals vs. Predicted Response',
                   y_label = 'PRESS Residuals')
norm_prob_plot(residual_var = PRESS_res, x_label = 'PRESS Residuals')
res_vs_regressor(residual_var = PRESS_res,
                 main_title1 = 'PRESS Residuals vs. Space time, min',
                 main_title2 = 'PRESS Residuals vs. Solvent total',
                 ylabel = 'PRESS Residuals')


### Problem 2
df <- MPV::table.b4
set.seed(2); chosen_rows <- sort(sample(seq(1,24), 15))
df <- df[chosen_rows,c(1, 5, 8, 10)]
n <- nrow(df)
ones <- rep(1, n)
y <- df[,1]
X <- cbind(ones, df[,c(2:4)])

beta_hat <- beta_hat_calc(X=X, y=y)
H <- H_calc(X)
y_hat <- y_hat_calc(H=H, y=y)
e <- e_calc(y=y, y_hat=y_hat)
### part (a)
SS_T_calc <- function(y) {
  n <- length(y)
  SS_T <- (t(y) %*% y) - ((sum(y)^2) / n)
  return(SS_T)
}
SS_R_calc <- function(beta_hat, X, y) {
  n <- length(y)
  SS_R <- (t(beta_hat) %*% t(X) %*% y) - ((sum(y)^2) / n)
  return(SS_R)
}

SS_Res <- SS_Res_calc(beta_hat = beta_hat, X = X, y = y)
SS_T <- SS_T_calc(y = y)
SS_R <- SS_R_calc(beta_hat = beta_hat, X = X, y = y)

SS_Res == SS_T - SS_R
SS_Res; SS_T;  SS_R
k <- 3; p <- k + 1
MS_R <- SS_R / k
MS_Res <- SS_Res / (n - k - 1)
F_0 <- MS_R / MS_Res

alpha <- 0.05
qf(p = (1 - alpha), df1 = k, df2 = (n - k - 1))
pf(q = F_0, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

r_squared_calc <- function(SS_R, SS_T) {
  r_squared <- SS_R / SS_T
  return(r_squared)
}
adj_r_squared_calc <- function(SS_Res, SS_T, n, p) {
  adj_r_squared <- 1 - ((SS_Res / (n - p)) / (SS_T / (n - 1)))
  return(adj_r_squared)
}

r_squared <- r_squared_calc(SS_R = SS_R, SS_T = SS_T)
adj_r_squared <- adj_r_squared_calc(
  SS_Res = SS_Res, SS_T = SS_T, n = n, p = p)

# ordinary residuals
par(mfrow = c(2,2))
norm_prob_plot(residual_var = e, x_label = 'Sorted Residuals')
res_vs_fitted_plot(residual_var = e,
                   main_title = 'Residuals vs. Predicted Response',
                   y_label = 'Residuals')

order(e, decreasing = TRUE)

# studentized residuals
H_diag <- diag(H)
studentized_residuals <- sapply(1:n, function(x) {
  e[x] / sqrt(MS_Res * (1 - H_diag[x]))
})
norm_prob_plot(residual_var = studentized_residuals, x_label = 'Studentized Residuals')
res_vs_fitted_plot(residual_var = studentized_residuals,
                   main_title = 'Studentized Residuals vs. Predicted Response',
                   y_label = 'Studentized Residuals')

# R-student residuals
par(mfrow = c(2,2))
S_squared <- sapply(1:n, function(x) {
  ((n - p) * MS_Res - ((e[x]^2) / (1 - H_diag[x]))) / (n - p - 1)
})

R_student_res <- sapply(1:n, function(x) {
  e[x] / sqrt(S_squared[x] * (1 - H_diag[x]))
})

norm_prob_plot(residual_var = R_student_res, x_label = 'R-student')
res_vs_fitted_plot(residual_var = R_student_res,
                   main_title = 'R-student vs. Predicted Response',
                   y_label = 'R-student')

# standardized residuals
standardized_res <- sapply(1:n, function(x) e[x] / sqrt(MS_Res))
norm_prob_plot(residual_var = standardized_res, x_label = 'Standardized Residuals')
res_vs_fitted_plot(residual_var = standardized_res,
                   main_title = 'Standardized Residuals vs. Predicted Response',
                   y_label = 'Standardized Residuals')

# PRESS residuals
par(mfrow = c(2,2))
PRESS_res <- sapply(1:n, function(x) e[x] / (1 - H_diag[x]))
norm_prob_plot(residual_var = PRESS_res, x_label = 'PRESS Residuals')
res_vs_fitted_plot(residual_var = PRESS_res,
                   main_title = 'PRESS Residuals vs. Predicted Response',
                   y_label = 'PRESS Residuals')

x4 <- X[,2]; x7 <- X[,3]; x9 <- X[,4]
par(mfrow = c(2,2))
plot(x4, e,
     main = 'Residuals vs. Living Space (sq ft x 1000)',
     ylab = 'Residuals',
     xlab = 'Living Space (sq ft x 1000)')
abline(h = 0)

plot(x7, e,
     main = 'Residuals vs. Number of Bedrooms',
     ylab = 'Residuals',
     xlab = 'Number of Bedrooms')
abline(h = 0)

plot(x9, e,
     main = 'Residuals vs. Number of Fireplaces',
     ylab = 'Residuals',
     xlab = 'Number of Fireplaces')
abline(h = 0)


### part (b)
data.frame(table(y))

### Problem 3
### part (a)
df <- MPV::p5.5
set.seed(3); chosen_rows <- sort(sample(seq(1,14), 7))
df <- df[chosen_rows,]
y <- df[,1]
n <- nrow(df)
ones <- rep(1, n)
X <- cbind(ones, df[,2])

beta_hat <- beta_hat_calc(X=X, y=y)
H <- H_calc(X)
y_hat <- y_hat_calc(H=H, y=y)
e <- e_calc(y=y, y_hat=y_hat)

SS_Res <- SS_Res_calc(beta_hat = beta_hat, X = X, y = y)
SS_T <- SS_T_calc(y = y)
SS_R <- SS_R_calc(beta_hat = beta_hat, X = X, y = y)

SS_Res == SS_T - SS_R
SS_Res; SS_T;  SS_R
k <- 1; p <- k + 1
MS_R <- SS_R / k
MS_Res <- SS_Res / (n - k - 1)
F_0 <- MS_R / MS_Res

alpha <- 0.01
qf(p = (1 - alpha), df1 = k, df2 = (n - k - 1))
pf(q = F_0, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

r_squared <- r_squared_calc(SS_R = SS_R, SS_T = SS_T)

# ordinary residuals
par(mfrow = c(2,2))
norm_prob_plot(residual_var = e, x_label = 'Sorted Residuals')
res_vs_fitted_plot(residual_var = e,
                   main_title = 'Residuals vs. Predicted Response',
                   y_label = 'Residuals')

order(e, decreasing = TRUE)

# studentized residuals
H_diag <- diag(H)
studentized_residuals <- sapply(1:n, function(x) {
  e[x] / sqrt(MS_Res * (1 - H_diag[x]))
})
norm_prob_plot(residual_var = studentized_residuals, x_label = 'Studentized Residuals')
res_vs_fitted_plot(residual_var = studentized_residuals,
                   main_title = 'Studentized Residuals vs. Predicted Response',
                   y_label = 'Studentized Residuals')

# R-student residuals
par(mfrow = c(2,2))
S_squared <- sapply(1:n, function(x) {
  ((n - p) * MS_Res - ((e[x]^2) / (1 - H_diag[x]))) / (n - p - 1)
})

R_student_res <- sapply(1:n, function(x) {
  e[x] / sqrt(S_squared[x] * (1 - H_diag[x]))
})

norm_prob_plot(residual_var = R_student_res, x_label = 'R-student')
res_vs_fitted_plot(residual_var = R_student_res,
                   main_title = 'R-student vs. Predicted Response',
                   y_label = 'R-student')

# standardized residuals
standardized_res <- sapply(1:n, function(x) e[x] / sqrt(MS_Res))
norm_prob_plot(residual_var = standardized_res, x_label = 'Standardized Residuals')
res_vs_fitted_plot(residual_var = standardized_res,
                   main_title = 'Standardized Residuals vs. Predicted Response',
                   y_label = 'Standardized Residuals')

# PRESS residuals
par(mfrow = c(2,2))
PRESS_res <- sapply(1:n, function(x) e[x] / (1 - H_diag[x]))
norm_prob_plot(residual_var = PRESS_res, x_label = 'PRESS Residuals')
res_vs_fitted_plot(residual_var = PRESS_res,
                   main_title = 'PRESS Residuals vs. Predicted Response',
                   y_label = 'PRESS Residuals')

par(mfrow = c(2,2))
plot(X[,2], e,
     main = 'Residuals vs. Weeks',
     ylab = 'Residuals',
     xlab = 'Weeks')
abline(h = 0)

plot(X[,2], studentized_residuals,
     main = 'Studentized Residuals vs. Weeks',
     ylab = 'Studentized Residuals',
     xlab = 'Weeks')
abline(h = 0)

plot(X[,2], R_student_res,
     main = 'R-Student vs. Weeks',
     ylab = 'R-Student',
     xlab = 'Weeks')
abline(h = 0)

plot(X[,2], standardized_res,
     main = 'Standardized Residuals vs. Weeks',
     ylab = 'Standardized Residuals',
     xlab = 'Weeks')
abline(h = 0)

par(mfrow = c(1,1))
plot(X[,2], PRESS_res,
     main = 'PRESS Residuals vs. Weeks',
     ylab = 'PRESS Residuals',
     xlab = 'Weeks')
abline(h = 0)

### part (b)
par(mfrow = c(2,2))
plot(df$weeks, df$defects,
     main = 'Defects per 10,000 vs. Weeks',
     xlab = 'Weeks', ylab = 'Defects per 10,000')
abline(beta_hat)

log_y <- log(df$defects)
log_x <- log(df$weeks)
plot(log_x, log_y,
     main = 'log Defects per 10,000 vs. log Weeks',
     xlab = 'log Weeks', ylab = 'log Defects per 10,000')
X_log <- cbind(ones, log_x)
beta_hat_log <- beta_hat_calc(X=X_log, y=log_y)
H_log <- H_calc(X_log)
y_hat_log <- y_hat_calc(H=H_log, y=log_y)
e_log <- e_calc(y=log_y, y_hat=y_hat_log)
abline(beta_hat_log)

SS_Res <- SS_Res_calc(beta_hat = beta_hat_log, X = X_log, y = log_y)
SS_T <- SS_T_calc(y = log_y)
SS_R <- SS_R_calc(beta_hat = beta_hat_log, X = X_log, y = log_y)

SS_Res == SS_T - SS_R
SS_Res; SS_T;  SS_R
k <- 1; p <- k + 1
MS_R <- SS_R / k
MS_Res <- SS_Res / (n - k - 1)
F_0 <- MS_R / MS_Res

alpha <- 0.01
qf(p = (1 - alpha), df1 = k, df2 = (n - k - 1))
pf(q = F_0, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

r_squared <- r_squared_calc(SS_R = SS_R, SS_T = SS_T)

# ordinary residuals
par(mfrow = c(2,2))
e <- e_log
norm_prob_plot(residual_var = e, x_label = 'Sorted Residuals')
res_vs_fitted_plot(residual_var = e,
                   main_title = 'Residuals vs. Predicted Response',
                   y_label = 'Residuals')

order(e, decreasing = TRUE)

# studentized residuals
H_diag <- diag(H)
studentized_residuals <- sapply(1:n, function(x) {
  e[x] / sqrt(MS_Res * (1 - H_diag[x]))
})
norm_prob_plot(residual_var = studentized_residuals, x_label = 'Studentized Residuals')
res_vs_fitted_plot(residual_var = studentized_residuals,
                   main_title = 'Studentized Residuals vs. Predicted Response',
                   y_label = 'Studentized Residuals')

# R-student residuals
par(mfrow = c(2,2))
S_squared <- sapply(1:n, function(x) {
  ((n - p) * MS_Res - ((e[x]^2) / (1 - H_diag[x]))) / (n - p - 1)
})

R_student_res <- sapply(1:n, function(x) {
  e[x] / sqrt(S_squared[x] * (1 - H_diag[x]))
})

norm_prob_plot(residual_var = R_student_res, x_label = 'R-student')
res_vs_fitted_plot(residual_var = R_student_res,
                   main_title = 'R-student vs. Predicted Response',
                   y_label = 'R-student')

# standardized residuals
standardized_res <- sapply(1:n, function(x) e[x] / sqrt(MS_Res))
norm_prob_plot(residual_var = standardized_res, x_label = 'Standardized Residuals')
res_vs_fitted_plot(residual_var = standardized_res,
                   main_title = 'Standardized Residuals vs. Predicted Response',
                   y_label = 'Standardized Residuals')

# PRESS residuals
par(mfrow = c(2,2))
PRESS_res <- sapply(1:n, function(x) e[x] / (1 - H_diag[x]))
norm_prob_plot(residual_var = PRESS_res, x_label = 'PRESS Residuals')
res_vs_fitted_plot(residual_var = PRESS_res,
                   main_title = 'PRESS Residuals vs. Predicted Response',
                   y_label = 'PRESS Residuals')

par(mfrow = c(2,2))
plot(X[,2], e,
     main = 'Residuals vs. log Weeks',
     ylab = 'Residuals',
     xlab = 'log Weeks')
abline(h = 0)

plot(X[,2], studentized_residuals,
     main = 'Studentized Residuals vs. log Weeks',
     ylab = 'Studentized Residuals',
     xlab = 'log Weeks')
abline(h = 0)

plot(X[,2], R_student_res,
     main = 'R-Student vs. log Weeks',
     ylab = 'R-Student',
     xlab = 'log Weeks')
abline(h = 0)

plot(X[,2], standardized_res,
     main = 'Standardized Residuals vs. log Weeks',
     ylab = 'Standardized Residuals',
     xlab = 'log Weeks')
abline(h = 0)

par(mfrow = c(1,1))
plot(X[,2], PRESS_res,
     main = 'PRESS Residuals vs. log Weeks',
     ylab = 'PRESS Residuals',
     xlab = 'log Weeks')
abline(h = 0)



