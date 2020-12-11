library(MPV)
### 5
### a
set.seed(1)
n <- 1e3
x1 <- rnorm(n = n, mean = 0, sd = 1)
x2 <- rnorm(n = n, mean = 3, sd = 0.2)
error <- rnorm(n = n, mean = 0, sd = 0.1)
y <- 0.2 * x1 + 3 * x2 + error

beta_hat_calc <- function(X, y) {
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
  return(beta_hat)
}
H_calc <- function(X) {
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

# full model
ones <- rep(1, n)
X_full <- cbind(ones, x1, x2)
beta_hat_full <- beta_hat_calc(X = X_full, y = y)
H_full <- H_calc(X = X_full)
y_hat_full <- y_hat_calc(H = H_full, y = y)
e_full <-  e_calc(y = y, y_hat = y_hat_full)
t(e_full) %*% e_full # 10.59006

# SLR model
ones <- rep(1, n)
X_red <- cbind(ones, x1)
beta_hat_red <- beta_hat_calc(X = X_red, y = y)
H_red <- H_calc(X = X_red)
y_hat_red <- y_hat_calc(H = H_red, y = y)
e_red <-  e_calc(y = y, y_hat = y_hat_red)
t(e_red) %*% e_red # 402.4067


### b

### 6
set.seed(1); n <- 22 # Sample rows
random_rows <- sort(sample(x = seq(1, 32), size = n, replace = FALSE))
na_rows <- c(23, 25)
random_rows <- random_rows[!random_rows %in% na_rows]
random_rows <- c(random_rows, 3, 6) # Ignore NA rows

# Create table subset
table_b3 <- MPV::table.b3; car_data <- table_b3[random_rows,]
car_data$index <- as.numeric(row.names(car_data))
car_data <- car_data[order(car_data$index),]
car_data$index <- NULL # Fix the data set

### a
# Subset data
y <- car_data[,1]; x1 <- car_data[,2]; x6 <- car_data[,7]; ones <- rep(1, n)
X <- cbind(ones, x1, x6)
beta_hat <- beta_hat_calc(X = X, y = y)
H <- H_calc(X = X)
y_hat <- y_hat_calc(H = H, y = y)
# e <- e_calc(y = y, y_hat = y_hat)
# t(e) %*% e # 148.2573

### b
SS_Res_calc <- function(y, beta_hat, X) {
  SS_Res <- (t(y) %*% y) - (t(beta_hat) %*% t(X) %*% y)
  return(SS_Res)
}
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
SS_Res; SS_T;  SS_R # 148.2573, 874.7109, 726.4536
k <- 2; n - k - 1; n - 1
MS_R <- SS_R / k # 363.2268
MS_Res <- SS_Res / (n - k - 1) # 7.803015
F <- MS_R / MS_Res # 46.54955

alpha <- 0.01
qf(p = (1 - alpha), df1 = k, df2 = (n - k - 1))
pf(q = F, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

### c
p <- k + 1
r_squared_calc <- function(SS_R, SS_T) {
  r_squared <- SS_R / SS_T
  return(r_squared)
}
adj_r_squared_calc <- function(SS_Res, SS_T, n, p) {
  adj_r_squared <- 1 - ((SS_Res / (n - p)) / (SS_T / (n - 1)))
  return(adj_r_squared)
}

r_squared <- r_squared_calc(SS_R = SS_R, SS_T = SS_T) # 0.8305071
adj_r_squared <- adj_r_squared_calc(
  SS_Res = SS_Res, SS_T = SS_T, n = n, p = p) # 0.8126658

### 2.4
k <- 1; p <- k + 1
X <- cbind(ones, x1)
beta_hat <- beta_hat_calc(X = X, y = y)
H <- H_calc(X = X)
y_hat <- y_hat_calc(H = H, y = y)
SS_Res <- SS_Res_calc(beta_hat = beta_hat, X = X, y = y)
SS_T <- SS_T_calc(y = y)
SS_R <- SS_R_calc(beta_hat = beta_hat, X = X, y = y)

r_squared <- r_squared_calc(SS_R = SS_R, SS_T = SS_T) # 0.7911184
adj_r_squared <- adj_r_squared_calc(
  SS_Res = SS_Res, SS_T = SS_T, n = n, p = p) # 0.7806743

### d
beta_hat_1 <- beta_hat[2]
C <- solve(t(X) %*% X); C_11 <- C[2,2]
sigma_hat_squared <- MS_Res
alpha <- 0.05
t_stat <- qt(p = (1 - alpha / 2), df = n - p)

CI <- c(beta_hat_1 - t_stat * sqrt(sigma_hat_squared * C_11),
        beta_hat_1 + t_stat * sqrt(sigma_hat_squared * C_11))

### e
test_stat_1 <- beta_hat_1 / sqrt(sigma_hat_squared * C_11) # -8.793749

beta_hat_6 <- beta_hat[3]; C_66 <- C[3,3]
test_stat_6 <- beta_hat_6 / sqrt(sigma_hat_squared * C_66) # 2.101295

t_stat <- qt(p = (1 - alpha / 2), df = (n - k - 1))

### f
x_0 <- c(1, 275, 2)
y_hat_0 <- t(x_0) %*% beta_hat
var_y_hat_0 <- sigma_hat_squared * (t(x_0) %*% C %*% x_0)
test_stat_0 <- qt(p = (1 - alpha / 2), df = n - p)

CI_0 <- c(y_hat_0 - test_stat_0 * sqrt(var_y_hat_0),
          y_hat_0 + test_stat_0 * sqrt(var_y_hat_0))

### g
x_0 <- c(1, 257, 2)
y_hat_0 <- t(x_0) %*% beta_hat
var_y_hat_0 <- sigma_hat_squared * (1 + (t(x_0) %*% C %*% x_0))
test_stat_0 <- qt(p = (1 - alpha / 2), df = n - p)

PI_0 <- c(y_hat_0 - test_stat_0 * sqrt(var_y_hat_0),
          y_hat_0 + test_stat_0 * sqrt(var_y_hat_0))
