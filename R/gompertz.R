#' Probability density function for the Gompertz distribution, with Median parametrization.
#'
#' @param x Value
#' @param mu Median parameter
#' @param beta Scale parameter
#' @param log Optional argument. If TRUE, returns log(pdf).
#'
#' @return PDF of gompertz distribution, with median parametrization.
#' @export
#'
#' @examples x <- seq(from = 0, to = 5, length.out = 100)
#' y <- bayesim::dgompertz(x, mu = 10, beta = 10)
#' plot(x, y, type = "l", ylab = "Density", main = "dgompertz(mu=10, beta=10)")
#' # Compare to online ressources
dgompertz <- function(x, mu, beta, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("dgompertz is only defined for x > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("dgompertz is only defined for mu > 0")
  }
  if (isTRUE(beta <= 0)) {
    stop("dgompertz is only defined for beta > 0")
  }
  lpdf <- log(-beta * log(0.5)) -
    log(exp(mu * beta) - 1) +
    (beta * x + (log(0.5) / (exp(mu * beta) - 1)) * (exp(beta * x) - 1))
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function for the Gompertz distribution, with Median parametrization.
#'
#' @param p Quantile to be calculated
#' @param mu Median parameter
#' @param beta Scale parameter
#'
#' @return Inverse of CDF, calculates a value, given a probability p
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, bayesim::qgompertz(x, mu = 10, beta = 1), type = "l", ylab = "Quantile", main = "apex-after-origin Gompertz(mu=10, beta=1)")
qgompertz <- function(p, mu, beta) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("qgompertz is only defined for 0 < p < 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("qgompertz is only defined for mu > 0")
  }
  if (isTRUE(beta <= 0)) {
    stop("qgompertz is only defined for beta > 0")
  }
  return(
    log1p(((exp(mu * beta) - 1) / (log(0.5))) * log1p(-p)) / beta
  )
}

#' RNG function for the Gompertz distribution, with Median parametrization.
#'
#' @param n Number of draws
#' @param mu Median parameter
#' @param beta Scale parameter
#'
#' @return A Gompertz distributed RNG vector of size n
#' @export
#'
#' @examples y <- bayesim::rgompertz(n, mu = 2, eta = 0.1)
#' hist(y, main = c(paste("Median:", median(y)), " for RNG of apex-after-origin Gompertz(mu=2, eta=0.1)"))
rgompertz <- function(n, mu, beta) {
  if (isTRUE(mu <= 0)) {
    stop("rgompertz is only defined for mu > 0")
  }
  if (isTRUE(beta <= 0)) {
    stop("rgompertz is only defined for beta > 0")
  }
  return(
    log1p(((exp(mu * beta) - 1) / (log(0.5))) * log(runif(n))) / beta
  )
}


#' Log-Likelihood vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of gompertz given data in prep
#'
#' @examples
log_lik_gompertz <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dgompertz(y, mu, beta, log = TRUE))
}

#' Posterior-Prediction vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of gompertz, given data in prep
#'
#' @examples
posterior_predict_gompertz <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  return(rgompertz(prep$ndraws, mu, beta))
}

#' Expectation-Predict vignette for the Gompertz distribution, with Median parametrization.
#' Not defined for the Gompertz family.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_gompertz <- function(prep) {
  stop("posterior_epred is not defined for the gompertz family")
}


#' Gompertz Stan-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_eta Link function for eta argument
#'
#' @return BRMS gompertz distribution family
#' @export
#'
#' @examples data <- list(a = a, y = bayesim::rgompertz(n, exp(0.5 * a + 1), 0.2))
#' fit1 <- brm(y ~ 1 + a,
#'   data = data, family = bayesim::gompertz(),
#'   stanvars = bayesim::gompertz()$stanvars, backend = "cmdstan"
#' )
#' plot(fit1)
gompertz <- function(link = "log", link_b = "log") {
  family <- brms::custom_family(
    "gompertz",
    dpars = c("mu", "beta"),
    links = c(link, link_b),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_gompertz,
    posterior_predict = posterior_predict_gompertz,
    posterior_epred = posterior_epred_gompertz
  )

  family$stanvars <- brms::stanvar(
    scode = "
      real gompertz_lpdf(real y, real mu, real beta) {
        return(
          log(-beta * log(0.5)) -
          log(exp(mu * beta) - 1) +
          (beta * y + (log(0.5) / (exp(mu * beta) - 1)) * (exp(beta * y) - 1))
        );
      }
      real gompertz_rng(real mu, real beta) {
      return(
        log1p(((exp(mu * beta) - 1) / (log(0.5))) * log(uniform_rng(0,1))) / beta
      );
      }",
    block = "functions"
  )
  return(family)
}
