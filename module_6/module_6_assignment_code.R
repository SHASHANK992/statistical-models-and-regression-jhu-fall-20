library(MPV)

### Problem 7.6
orig_df <- MPV::p7.6
n <- nrow(orig_df)
ones <- rep(1, n)
y <- orig_df[,1]
x1 <- orig_df[,2]
x2 <- orig_df[,3]

# part (a) Fit a second-order polynomial.
df_a <- cbind(ones, x1, x2, x1^2, x2^2, x1*x2)

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

beta_hat_a <- beta_hat_calc(X = df_a, y = y)
H_a <- H_calc(X = df_a)
y_hat_a <- y_hat_calc(H = H_a, y = y)
e_a <-  e_calc(y = y, y_hat = y_hat_a)
t(e_a) %*% e_a # 2.302162

# part (b) Test for significance of regression.
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

SS_Res_a <- SS_Res_calc(beta_hat = beta_hat_a, X = df_a, y = y)
SS_T_a <- SS_T_calc(y = y)
SS_R_a <- SS_R_calc(beta_hat = beta_hat_a, X = df_a, y = y)

SS_Res_a == SS_T_a - SS_R_a
SS_Res_a; SS_T_a;  SS_R_a # 2.302142, 342.1899, 339.8878
k <- 5; n - k - 1; n - 1
MS_R_a <- SS_R_a / k # 67.97755
MS_Res_a <- SS_Res_a / (n - k - 1) # 0.3836904
F_a <- MS_R_a / MS_Res_a # 177.1677

alpha <- 0.1
qf(p = (1 - alpha), df1 = k, df2 = (n - k - 1))
pf(q = F_a, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

# Double-check with lm / anova
m1 <- lm(y ~ x1*x2 + I(x1^2) + I(x2^2))
m2 <- lm(y~1)
anova(m1, m2)

# part (d) Does the interaction term contribute significantly to the model?
C_mat_calc <- function(X) {
  C_mat <- solve(t(X) %*% X)
  return(C_mat)
}

sigma_hat_squared <- MS_Res_a
C_mat_a <- C_mat_calc(X = df_a)
t_a <- (beta_hat_a[6] / sqrt(sigma_hat_squared * C_mat_a[6,6]))

alpha <- 0.05
qt(p = (1 - alpha / 2), df = (n - k - 1))

# Reference: https://stats.stackexchange.com/questions/45153/manually-calculating-p-value-from-t-value-in-t-test
2 * pt(q = abs(t_a), df = (n - k - 1), lower.tail = FALSE)

# Double-check with summary
summary(m1)

# part (e) Do the second-order terms contribute significantly to the model?
df_a1 <- df_a[,c(1,2,3)]
beta_hat_red <- beta_hat_calc(X = df_a1, y = y)
SS_R_b1_a <- SS_R_calc(beta_hat = beta_hat_red, X = df_a1, y = y)
SS_R_b2_given_b1_a <- SS_R_a - SS_R_b1_a

r <- 3
F_a1 <- (SS_R_b2_given_b1_a / r) / MS_Res_a # 5.057906
alpha <- 0.05
qf(p = (1 - alpha), df1 = r, df2 = (n - k - 1))
pf(q = F_a1, df1 = r, df2 = (n - k - 1), lower.tail = FALSE)

# Double-check
m2 <- lm(y ~ x1 + x2)
anova(m1, m2)

### Problem 8.11
orig_df <- MPV::p8.11
n <- nrow(orig_df)
y <- orig_df[,1]
Xs <- orig_df[,2]
ones <- rep(1, n)

indicator_ones <- rep(1, 5)
indicator_zeros <- rep(0, 5)
x1 <- c(indicator_ones, rep(indicator_zeros, 4))
x2 <- c(rep(indicator_zeros, 1), indicator_ones, rep(indicator_zeros, 3))
x3 <- c(rep(indicator_zeros, 2), indicator_ones, rep(indicator_zeros, 2))
x4 <- c(rep(indicator_zeros, 3), indicator_ones, rep(indicator_zeros, 1))

X <- cbind(ones, x1, x2, x3, x4)

# part (b)	Find the least-squares estimates of the model parameters.
beta_hat <- beta_hat_calc(X = X, y = y)
H <- H_calc(X = X)
y_hat <- y_hat_calc(H = H, y = y)
e <-  e_calc(y = y, y_hat = y_hat)
t(e) %*% e # 161.2

# Double-check
m1 <- lm(y~X-1)

# part (c)	Find a point estimate of the difference in mean strength between 15% and 25% cotton.
beta_hat[2] - beta_hat[4] # -7.8

# part (d)	Test the hypothesis that the mean tensile strength is the same for all five cotton percentages.
k <- 4
SS_R_full <- SS_R_calc(beta_hat = beta_hat, X = X, y = y)
MS_R_full <- SS_R_full / k # 118.94
SS_Res <- SS_Res_calc(beta_hat = beta_hat, X = X, y = y)
MS_Res <- SS_Res / (n - k - 1) # 8.06
F_0 <- MS_R_full / MS_Res # 14.75682
pf(q = F_0, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

# Double-check
X_red <- X[,1]
m2 <- lm(y~X_red-1)
anova(m1, m2)

### Problem 8.16
Location <- seq(1, 17)
INHIBIT <- c(0.00, 1.00, 6.00, 7.00, 7.00, 7.00, 9.00, 9.50, 10.00, 11.00, 12.50, 14.00, 20.00, 21.00, 25.00, 39.00, 59.00)
UVB <- c(0.00, 0.00, 0.01, 0.01, 0.02, 0.03, 0.04, 0.01, 0.00, 0.03, 0.03, 0.01, 0.03, 0.04, 0.02, 0.03, 0.03)
SURFACE <- c('Deep', 'Deep', 'Deep', 'Surface', 'Surface', 'Surface', 'Surface', 'Deep', 'Deep', 'Surface', 'Surface', 'Deep', 'Deep', 'Surface', 'Deep', 'Deep', 'Deep')
orig_df <- data.frame(Location, INHIBIT, UVB, SURFACE)

set.seed(1); chosen_rows <- sample(Location, 12)
chosen_rows <- sort(chosen_rows)
n <- length(chosen_rows)
ones <- rep(1, n)
df <- orig_df[chosen_rows,]
y <- df$INHIBIT
x1 <- df$UVB

# Reference: https://stackoverflow.com/questions/40780088/r-code-categorical-variable-to-1-and-0
df$SURFACE <- as.factor(df$SURFACE)
df$is_Deep <- as.numeric(df$SURFACE)
df[df$is_Deep == 2,]$is_Deep <- 0
x2 <- df$is_Deep
X <- cbind(ones, x1, x2, x1*x2)

# Test for significance of regression
k <- 3
beta_hat <- beta_hat_calc(X = X, y = y)
SS_R_full <- SS_R_calc(beta_hat = beta_hat, X = X, y = y)
MS_R_full <- SS_R_full / k
SS_Res <- SS_Res_calc(beta_hat = beta_hat, X = X, y = y)
MS_Res <- SS_Res / (n - k - 1)
F_0 <- MS_R_full / MS_Res # 9.119374
pf(q = F_0, df1 = k, df2 = (n - k - 1), lower.tail = FALSE)

# Double-check
m1 <- lm(y ~ x1*x2); m2 <- lm(y ~ 1)
anova(m1, m2)

# Test for interaction term significance
sigma_hat_squared <- MS_Res
C_mat_a <- C_mat_calc(X = X)
t_0 <- (beta_hat[4] / sqrt(sigma_hat_squared * C_mat_a[4,4]))

alpha <- 0.05
qt(p = (1 - alpha / 2), df = (n - k - 1))

2 * pt(q = abs(t_0), df = (n - k - 1), lower.tail = FALSE)

# Double-check with summary
summary(m1)

# Test b1 and b2
t_0_b1 <- (beta_hat[2] / sqrt(sigma_hat_squared * C_mat_a[2,2]))
t_0_b2 <- (beta_hat[3] / sqrt(sigma_hat_squared * C_mat_a[3,3]))

2 * pt(q = abs(t_0_b1), df = (n - k - 1), lower.tail = FALSE)
2 * pt(q = abs(t_0_b2), df = (n - k - 1), lower.tail = FALSE)

