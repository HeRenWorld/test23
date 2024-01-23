rm(list = ls())
cat("\014")
Sys.setenv(LANG = "en")

test test

# Prob.1
# 1.d
N <- 12
time <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
data.smk <- c(29, 16, 17, 4, 3, 9, 4, 5, 1, 1, 1, 3, 7)
data.nonsmk <- c(198, 107, 55, 38, 18, 22, 7, 9, 5, 3, 6, 6, 12)

# mle
mle.p <-
  function(data, time, n) {
    # create a function with the name my_function
    z <- qnorm(0.95)
    
    phat <-
      sum(data[1:n]) / (sum((time * data)[1:n]) + n * data[n + 1])
    fisher <-
      sum(data[1:n]) / phat ^ 2 + sum(data * (time - 1)) / ((1 - phat) ^ 2)
    var <- 1 / fisher
    ci <- c(phat - z * sqrt(var), phat + z * sqrt(var))
    
    res <- list("phat" = phat,
                "var" = var,
                "CI" = ci)
    return(res)
  }

# phat.smk <- sum(data.smk[1:N]) / (sum((T * data.smk)[1:12]) + N * data.smk[13])
# phat.nonsmk <- sum(data.nonsmk[1:N]) / (sum((T * data.nonsmk)[1:12]) + N * data.nonsmk[13])
# variance
mle.smk <- mle.p(data.smk, time, N)
mle.nonsmk <- mle.p(data.nonsmk, time, N)

# 1.f
prior_fun <- function(parameter) {
  a = parameter[1]
  b = parameter[2]
  res = (pbeta(0.15, a, b) - 0.05) ^ 2 + (pbeta(0.4, a, b) - 0.95) ^ 2
  return(res)
}
prior.phat <- optim(par = c(1, 1), fn = prior_fun)
alpha.hat = prior.phat$par[1]
beta.hat = prior.phat$par[2]
alpha1.hat = alpha.hat + sum(data.smk[1:N])
beta1.hat = beta.hat + sum((time - 1) * data.smk)
alpha2.hat = alpha.hat + sum(data.nonsmk[1:N])
beta2.hat = beta.hat + sum((time - 1) * data.nonsmk)


# posterior mean
pm.1 = alpha1.hat / (alpha1.hat + beta1.hat)
pm.2 = alpha2.hat / (alpha2.hat + beta2.hat)
ic.post1 = c(qbeta(0.05, alpha1.hat, beta1.hat),
             qbeta(0.95, alpha1.hat, beta1.hat))
ic.post2 = c(qbeta(0.05, alpha2.hat, beta2.hat),
             qbeta(0.95, alpha2.hat, beta2.hat))

# plot
x = seq(0, 1, by = 0.001)
pdf.p1 = dbeta(x, alpha1.hat, beta1.hat)
pdf.p2 = dbeta(x, alpha2.hat, beta2.hat)
plot(x, pdf.p1, ylim = c(0, max(pdf.p1, pdf.p2)), col = "blue")
lines(x, pdf.p2, col = "red")
legend(
  "topright",
  legend = c("Smoker", "NonSmoker"),
  fill = c("blue", "red")
)



# Prob.2
# 2.d
log.lik <- function(parameter, time, y, n) {
  a = parameter[1]
  b = parameter[2]
  lbeta.tb = lbeta(a, time + b)[1:n]
  l = sum(y[1:n] * (log(a / (b + time - 1))[1:n] + lbeta.tb - lbeta(a, b))) +
    y[n + 1] * (lbeta(a, b + n) - lbeta(a, b))
  return(-l)
}
mle.ab1 <- optim(
  par = c(1, 1),
  fn = log.lik,
  time = time,
  y = data.smk,
  n = N
)
mle.ab2 <- optim(
  par = c(1, 1),
  fn = log.lik,
  time = time,
  y = data.nonsmk,
  n = N
)

# Prob.3
# 3.b
bin.p <- function(data, time, n) {
  z <- qnorm(0.95)
  Yt <- c(data[1:n], 0)
  total <- c((time * data)[1:n], N * data[n + 1])
  glm.p <-
    glm(cbind(Yt, total - Yt) ~ 1, family = binomial(link = "logit"))
  
  phat <- 1 / (1 + exp(-glm.p$coefficients))
  fisher <-
    sum(data[1:n]) / phat ^ 2 + sum(data * (time - 1)) / ((1 - phat) ^ 2)
  var <- 1 / fisher
  ci.asy <- c(phat - z * sqrt(var), phat + z * sqrt(var))
  
  ci <- 1 / (1 + exp(-confint(glm.p)))
  res <- list(
    "model" = glm.p,
    "phat" = phat,
    "CI.asy" = ci.asy,
    "CI" = ci
  )
  return(res)
}

bin.p1 <- bin.p(data.smk, time, N)
bin.p2 <- bin.p(data.nonsmk, time, N)

# 3.c
set.seed(123)
B = 10000
p1.sample = rbeta(B, alpha1.hat, beta1.hat)
p2.sample = rbeta(B, alpha2.hat, beta2.hat)
diff.p12 = p1.sample - p2.sample
hist(diff.p12)

p.diffp12 = sum(diff.p12 < 0) / B

# 3.f
prob.y <- function(alpha.hat, beta.hat, time, n) {
  yt = lbeta(alpha.hat + 1, time + beta.hat - 1) - lbeta(alpha.hat, beta.hat)
  yn = lbeta(alpha.hat, n + beta.hat) - lbeta(alpha.hat, beta.hat)
  return(exp(c(yt[1:n], yn)))
}
prob.y1 = prob.y(alpha1.hat, beta1.hat, time, N)
prob.y2 = prob.y(alpha2.hat, beta2.hat, time, N)
sample1 = rmultinom(B, N, prob.y1)
sample2 = rmultinom(B, N, prob.y2)
hist(sample1)
hist(sample2)
