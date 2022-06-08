library(bayesim)
library(extraDistr)
library(brms)

get_a <- function(mu, b) {
  # eta = a / b
  a <- b * (-log(0.5) / (exp(mu * b) - 1))
  return(a)
}

n <- 10000   # number of testvalues

x = seq(from = 0 , to = 10 , length.out=n)

# PDF test
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dgompertz_b(x, mu = 2, b = 0.1), type="l", ylab = "Density", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::dgompertz_b(x, mu = 1, b = 1), type="l", ylab = "Density", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::dgompertz_b(x, mu = 2, b = 3), type="l", ylab = "Density", main = "no-apex Gompertz(mu=2, eta=3)")

#stop("PDF")

# PDF comparison to a reference
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::dgompertz_b(x, mu = 2, b = 0.1) - extraDistr::dgompertz(x, get_a(2, 0.1), 0.1), ylab = "Density difference", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::dgompertz_b(x, mu = 1, b = 1) - extraDistr::dgompertz(x, get_a(1, 1), 1), ylab = "Density difference", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::dgompertz_b(x, mu = 2, b = 3) - extraDistr::dgompertz(x, get_a(2, 3), 3), ylab = "Density difference", main = "no-apex Gompertz(mu=2, eta=3)")

#stop("PDF-diff")

# Quantile test
x = seq(from = 0 , to = 1 , length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qgompertz_b(x, mu = 2, b = 0.1), type="l", ylab = "Quantile", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::qgompertz_b(x, mu = 1, b = 1), type="l", ylab = "Quantile", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::qgompertz_b(x, mu = 2, b = 3), type="l", ylab = "Quantile", main = "no-apex Gompertz(mu=2, eta=3)")

#stop("Quantiles")

# Quantile comparison to a reference
x = seq(from = 0 , to = 1 , length.out=n)
layout(matrix(1:3, ncol = 3))
plot(x, bayesim::qgompertz_b(x, mu = 2, b = 0.1) - extraDistr::qgompertz(x, get_a(2, 0.1), 0.1), ylab = "Quantile difference", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
plot(x, bayesim::qgompertz_b(x, mu = 1, b = 1) - extraDistr::qgompertz(x, get_a(1, 1), 1), ylab = "Quantile difference", main = "apex-at-origin Gompertz(mu=1, eta=1)")
plot(x, bayesim::qgompertz_b(x, mu = 2, b = 3) - extraDistr::qgompertz(x, get_a(2, 3), 3), ylab = "Quantile difference", main = "no-apex Gompertz(mu=2, eta=3)")

#stop("Quantile-diff")

layout(matrix(1:3, ncol = 3))
y = bayesim::rgompertz_b(n, mu = 2, b = 0.1)
hist(y, main = c(paste("Median:", median(y)), " for RNG of apex-after-origin Gompertz(mu=2, eta=0.1)"))
y = bayesim::rgompertz_b(n, mu = 1, b = 1)
hist(y, main = c(paste("Median:", median(y)), " for RNG of apex-at-origin Gompertz(mu=1, eta=1)"))
y = bayesim::rgompertz_b(n, mu = 2, b = 3)
hist(y, main = c(paste("Median:", median(y)), " for RNG of no-apex Gompertz(mu=2, eta=3)"))

#stop("RNG")

# Now all the R code is tested, now test the BRMS family
n = 1000
a = rnorm(n)
data = list(a = a, y = bayesim::rgompertz_b(n, exp(0.5 * a + 1), b = 2))
layout(1)
hist(log(data$y))

fit1 <- brm(
  y ~ 1 + a,
  init = 0.1,
  data = data,
  family = bayesim::gompertz_b(),
  stanvars = bayesim::gompertz_b()$stanvars,
  backend = "cmdstan",
  cores = 4
)

summary(fit1)

plot(fit1)

brms::pp_check(fit1)

brms::conditional_effects(fit1)

