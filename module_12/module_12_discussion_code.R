logit <- function(p) {
  return(log(p / (1 - p)))
}

probit <- function(p) {
  return(qnorm(p = p))
}

log_log <- function(p) {
  return(log(-log(1 - p)))
}

p_seq <- seq(0, 1, length.out = 1e3)
plot(p_seq, logit(p_seq), type = 'l',
     main = 'Link Functions vs. p',
     xlab = 'p', ylab = 'Link Function', lty = 2, col = 'green')
lines(p_seq, probit(p_seq), lty = 2, col = 'red')
lines(p_seq, log_log(p_seq), lty = 2, col = 'blue')
legend("topleft", legend = c('Logit', 'Probit', 'Log-log'),
       col = c('green', 'red', 'blue'), lty = rep(2,3))
abline(h = 0, lty = 2)
abline(v = 0.5, lty = 2)
