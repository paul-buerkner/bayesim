# implementation for b instead of eta
#' Title
#'
#' @param x
#' @param mu
#' @param b
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dgompertz_b <- function(x, mu, b, log = FALSE) {
  eta <- -log(0.5) / (exp(mu * b) - 1)
  lpdf <- log(b) + log(eta) + (eta + b * x - eta * exp(b * x))

  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Title
#'
#' @param p
#' @param mu
#' @param b
#'
#' @return
#' @export
#'
#' @examples
qgompertz_b <- function(p, mu, b) {
  eta <- -log(0.5) / (exp(mu * b) - 1)
  x <- (1 / b) * log1p(-(1 / eta) * log1p(-p))
  return(x)
}

#' Title
#'
#' @param n
#' @param mu
#' @param b
#'
#' @return
#' @export
#'
#' @examples
rgompertz_b <- function(n, mu, b) {
  return(qgompertz_b(runif(n), mu, b))
}


#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#' @export
#'
#' @examples
log_lik_gompertz_b <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  b <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dgompertz_b(y, mu, b, log = TRUE))
}

#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
posterior_predict_gompertz_b <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  b <- brms::get_dpar(prep, "beta", i = i)
  return(rgompertz(prep$ndraws, mu, b))
}

#' Title
#'
#' @param link
#' @param link_b
#'
#' @return
#' @export
#'
#' @examples
gompertz_b <- function(link = "log", link_b = "log") {
  family <- brms::custom_family(
    "gompertz_b",
    dpars = c("mu", "beta"),
    links = c(link, link_b),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_gompertz_b,
    posterior_predict = posterior_predict_gompertz_b,
    posterior_epred = posterior_epred_gompertz
  )

  family$stanvars <- brms::stanvar(
    scode = "
      real gompertz_b_lpdf(real y, real mu, real beta) {

        real eta = -log(0.5) / (exp(mu * beta) - 1);
        real lpdf = log(beta) + log(eta) + (eta + beta * y - eta * exp(beta * y));
        return(lpdf);
      }
      real gompertz_b_rng(real mu, real beta) {

        real eta = -log(0.5) / (exp(mu * beta) - 1);
        real x = (1 / beta) * log1p(-(1 / eta) * log1p(-uniform_rng(0, 1)));
        return(x);
      }",
    block = "functions"
  )
  return(family)
}

# implementation with a, b pdf instead
#' Title
#'
#' @param x
#' @param mu
#' @param b
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dgompertz_ab <- function(x, mu, b, log = FALSE) {
  a <- -(b * log(0.5)) / (exp(mu * b) - 1)
  lpdf <- log(a) + b * x - (a / b) * (exp(b * x) - 1)

  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Title
#'
#' @param p
#' @param mu
#' @param b
#'
#' @return
#' @export
#'
#' @examples
qgompertz_ab <- function(p, mu, b) {
  a <- -(b * log(0.5)) / (exp(mu * b) - 1)
  x <- (1 / b) * log1p(-(b/a) * log1p(-p))
  return(x)
}

#' Title
#'
#' @param n
#' @param mu
#' @param b
#'
#' @return
#' @export
#'
#' @examples
rgompertz_ab <- function(n, mu, b) {
  return(qgompertz_ab(runif(n), mu, b))
}

#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#' @export
#'
#' @examples
log_lik_gompertz_ab <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  b <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dgompertz_ab(y, mu, b, log = TRUE))
}

#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
posterior_predict_gompertz_ab <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  b <- brms::get_dpar(prep, "beta", i = i)
  return(rgompertz_ab(prep$ndraws, mu, b))
}

#' Title
#'
#' @param link
#' @param link_b
#'
#' @return
#' @export
#'
#' @examples
gompertz_ab <- function(link = "log", link_b = "log") {
  family <- brms::custom_family(
    "gompertz_ab",
    dpars = c("mu", "beta"),
    links = c(link, link_b),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_gompertz_ab,
    posterior_predict = posterior_predict_gompertz_ab,
    posterior_epred = posterior_epred_gompertz
  )

  family$stanvars <- brms::stanvar(
    scode = "
      real gompertz_ab_lpdf(real y, real mu, real beta) {

        real a = -(beta * log(0.5)) / (exp(mu * beta) - 1);
        real lpdf = log(a) + beta * y - (a / beta) * (exp(beta * y) - 1);

        return(lpdf);
      }
      real gompertz_ab_rng(real mu, real beta) {

        real a = -(beta * log(0.5)) / (exp(mu * beta) - 1);
        real x = (1 / beta) * log1p(-(beta/a) * log1p(uniform_rng(-1, 0)));
        return(x);
      }",
    block = "functions"
  )
  return(family)
}
