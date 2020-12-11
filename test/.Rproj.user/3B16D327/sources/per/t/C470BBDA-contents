# Generate data
N <- 1e3
# Reference: https://stats.stackexchange.com/questions/70855/generating-random-variables-from-a-mixture-of-normal-distributions
set.seed(1); u_sample <- runif(n = N)
set.seed(1)
x1 <- sapply(1:N, function(x) {
  if (u_sample[x] > 0.6) {
    return(rnorm(n = 1, mean = 0, sd = 1))
  } else {
    return(rnorm(n = 1, mean = 2.5, sd = 2))
  }
})
set.seed(1)
x2 <- sapply(1:N, function(x) {
  if (u_sample[x] > 0.1) {
    return(rnorm(n = 1, mean = 1, sd = 2))
  } else {
    return(rnorm(n = 1, mean = 3, sd = 5))
  }
})
x1x2 <- x1*x2

y <- 1.5 + 2 * x1 + 3 * x2 + 0.5 * x1x2 + rnorm(n = N, mean = 0, sd = 1)
par(mfrow=c(2,2))
plot(density(y), main = 'Density Plot of y')
plot(density(x1), main = 'Density Plot of x1')
plot(density(x2), main = 'Density Plot of x2')
plot(density(x1x2), main = 'Density Plot of x1x2')

# full model
df <- as.data.frame(cbind(y, x1, x2, x1x2))
full_model <- lm(formula = y~., data = df)
intercept_only <- lm(y~1, data = df)
summary(full_model)

# all possible regressors
best <- regsubsets(x = y~., data = df, nvmax = 3)
res.sum <- summary(best)
p.m <- 2:4
aic <- N * log(res.sum$rss / N) + 2 * p.m
MS_Res <- res.sum$rss / (N - p.m)
data.frame(
  R2 = which.max(res.sum$rsq),
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic),
  AIC = which.min(aic),
  MSRes = which.min(MS_Res)
)

data.frame(
  R2 = round(res.sum$rsq[3],4),
  Adj.R2 = round(res.sum$adjr2[3],4),
  CP = round(res.sum$cp[3],4),
  BIC = round(res.sum$bic[3],4),
  AIC = round(aic[3],4),
  MSRes = round(MS_Res[3],4)
)

# forward selection
forward <- MASS::stepAIC(intercept_only,
                         scope = list(upper=full_model, lower=intercept_only),
                         direction = c('forward'))

# backward elimination
backward <- MASS::stepAIC(full_model,
                           # scope = list(upper=full_model, lower=intercept_only),
                           direction = c('backward'))

# stepwise regression
stepwise <- MASS::stepAIC(intercept_only,
                       scope = list(upper=full_model, lower=intercept_only),
                       direction = c('both'))




