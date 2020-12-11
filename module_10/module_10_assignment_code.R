### Problem 1
df <- MPV::table.b2
n <- 20
set.seed(1); chosen_rows <- sort(sample(seq(1, nrow(df)), n))
df1 <- df[chosen_rows,]

# Reference: https://stats.stackexchange.com/questions/347652/default-stepaic-in-r
# (a) forward selection
model_0 <- lm(y~1, data = df1)
model_1 <- lm(y~., data = df1)
forward1 <- MASS::stepAIC(model_0,
                          scope = list(upper=model_1, lower=model_0),
                          direction = c('forward'))

# (b) backward elimination
backward1 <- MASS::stepAIC(model_1,
                           # scope = list(upper=model_1, lower=model_0),
                           direction = c('backward'))

# (c) stepwise regression
step1 <- MASS::stepAIC(model_0,
                       scope = list(upper=model_1, lower=model_0),
                       direction = c('both'))

# (d) all possible regressions
best1 <- leaps::regsubsets(x = y~., data = df1, nvmax = 5)
best1_sum <- summary(best1)
p.m <- 2:6
MS_Res <- best1_sum$rss / (n-p.m)
data.frame(
  rsq = which.max(best1_sum$rsq),
  CP = which.min(best1_sum$cp),
  MSRes = which.min(MS_Res)
)
which.min(best1_sum$bic)
aic <- n * log(best1_sum$rss / n) + 2 * p.m
which.min(aic)

### Problem 2
df <- MPV::table.b1
n <- 20
set.seed(1); chosen_rows <- sort(sample(seq(1, nrow(df)), n))
df2 <- df[chosen_rows,]

# (a)
# PRESS residuals
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
X <- df2[,2:ncol(df2)]
ones <- rep(1, n)
X <- cbind(ones, X)
y <- df2$y
beta_hat <- beta_hat_calc(X=X, y=y)
H <- H_calc(X)
y_hat <- y_hat_calc(H=H, y=y)
e <- e_calc(y=y, y_hat=y_hat)
H_diag <- diag(H)

PRESS_res <- sapply(1:n, function(x) e[x] / (1 - H_diag[x]))
PRESS <- sum(PRESS_res^2)
sum(e^2)

# 2.1.1
data.frame(e=e,PRESS=PRESS_res)

# (b)
set.seed(1)
chosen_row_subset <- sort(sample(chosen_rows, size = length(chosen_rows)/2))
df2b <- df[chosen_row_subset,]
n <- 10
X <- df2b[,2:ncol(df2)]
ones <- rep(1, n)
X <- cbind(ones, X)
y <- df2b$y
beta_hat <- beta_hat_calc(X=X, y=y)
H <- H_calc(X)
y_hat <- y_hat_calc(H=H, y=y)
e <- e_calc(y=y, y_hat=y_hat)
H_diag <- diag(H)

PRESS_res_b <- sapply(1:n, function(x) e[x] / (1 - H_diag[x]))
PRESS_b <- sum(PRESS_res_b^2)
sum(e^2)

# 2.2.1
data.frame(e=e,PRESS=PRESS_res_b)

# predictive powers
deleted_rows <- chosen_rows[!(chosen_rows %in% chosen_row_subset)]
df2c <- df[deleted_rows,]
X2 <- df2c[,2:ncol(df2c)]
X2 <- cbind(ones, X2)
y_hat2 <- as.matrix(X2) %*% as.matrix(beta_hat)
y2 <- df2c$y
e2 <- e_calc(y=y2, y_hat=y_hat2)
PRESS_res_c <- sapply(1:n, function(x) e2[x] / (1 - H_diag2[x]))
PRESS_c <- sum(PRESS_res_c^2)
sum(e2^2)

