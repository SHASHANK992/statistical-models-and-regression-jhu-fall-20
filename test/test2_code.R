### Problem 3
df <- data.frame(
  y=c(7,8,5,4,2,10,9,10,8,8),
  x1=c(9,6,10,8,5,7,6,5,5,4),
  x2=rep(c(1,-1), each=5)
)
y <- df[,1]
n <- nrow(df)
ones <- rep(1, n)
x1 <- df[,2]; x2 <- df[,3]
X <- cbind(ones, x1, x2, x1*x2)

### part (a)
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

X_red <- cbind(ones, x1)
beta_hat_full <- beta_hat_calc(X = X, y = y)
beta_hat_red <- beta_hat_calc(X = X_red, y = y)

SS_R_calc <- function(beta_hat, X, y) {
  n <- length(y)
  SS_R <- (t(beta_hat) %*% t(X) %*% y) - ((sum(y)^2) / n)
  return(SS_R)
}
SS_Res_calc <- function(y, beta_hat, X) {
  SS_Res <- (t(y) %*% y) - (t(beta_hat) %*% t(X) %*% y)
  return(SS_Res)
}

SS_Res_full <- SS_Res_calc(y = y, beta_hat = beta_hat, X = X)
k <- 3; p <- k + 1
r <- 2
MS_Res <- SS_Res_full / (n - p)
SS_R_full <- SS_R_calc(beta_hat = beta_hat, X = X, y = y)
SS_R_red <- SS_R_calc(beta_hat = beta_hat_red, X = X_red, y = y)

F_0 <- ((SS_R_full - SS_R_red) / r) / MS_Res
alpha <- 0.05
qf(p = (1 - alpha), df1 = k, df2 = (n - k - 1))
pf(q = F_0, df1 = r, df2 = (n - k - 1), lower.tail = FALSE)

### part (b)
2 * beta_hat[3] + 10 * beta_hat[4] # -4.38551

z_value <- qnorm(alpha/2, lower.tail = FALSE)
x_01 <- c(1, 5, 1, 5)
x_02 <- c(1, 5, -1, -5)
y_hat_01 <- t(x_01) %*% beta_hat
y_hat_02 <- t(x_02) %*% beta_hat
sigma_squared <- 2

CI_bound <- z_value * sqrt(sigma_squared * t(x_01 - x_02) %*%
                             solve(t(X) %*% X) %*%
                             (x_01 - x_02))

(y_hat_01 - y_hat_02) + CI_bound
(y_hat_01 - y_hat_02) - CI_bound

### part (c)
PI_bound <- z_value * sqrt(sigma_squared *
                             (2 +
                                t(x_01 - x_02) %*%
                                solve(t(X) %*% X) %*%
                                (x_01 - x_02)))

(y_hat_01 - y_hat_02) + PI_bound
(y_hat_01 - y_hat_02) - PI_bound

### part (d)
X_d <- X[,c(1,3)]
beta_hat_d <- beta_hat_calc(X = X_d, y = y)
H_d <- H_calc(X = X_d)
y_hat_d <- y_hat_calc(H = H_d, y = y)
e_d <- e_calc(y = y, y_hat = y_hat_d)

e_d[8]

residual_variance <- sigma_squared * (diag(n) - H_d)
residual_variance[8,8]
