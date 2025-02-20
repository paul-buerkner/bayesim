#' Title
#'
#' @param family_list
#'
#' @return
#' @export
#'
#' @examples
get_prefit <- function(family_list, brms_backend) {
  family <- brms_family_lookup(
    family_list$fit_family,
    family_list$fit_link
  )
  prefit <- brms::brm(
    y ~ 1 + x,
    data = list(y = c(0.5), x = c(1)),
    family = family,
    stanvars = family$stanvars,
    chains = 0,
    refresh = 0,
    silent = 2,
    backend = brms_backend,
    prior = prior_lookup(family_list$fit_family),
    init = 0.1
  )
  return(prefit)
}


#' Title
#'
#' @param fit_configuration
#'
#' @return
#' @export
#'
#' @examples
build_prefit_list <- function(fit_configuration, brms_backend) {
  prefit_configurations <- unique(fit_configuration[c("fit_family", "fit_link")])
  prefit_configurations <- lapply(split(
    prefit_configurations,
    sort(as.numeric(rownames(prefit_configurations)))
  ), as.list)

  results <- vector(mode = "list", length = length(prefit_configurations))
  for (i in seq_along(prefit_configurations)) {
    results[[i]] <- get_prefit(family_list = prefit_configurations[[i]], brms_backend = brms_backend)
  }

  prefit_list <- list()
  for (i in seq_along(results)) {
    prefit <- results[[i]]
    if (!is.null(prefit$family$name)) {
      family <- prefit$family$name
    } else {
      family <- prefit$family$family
    }
    key <- paste0(family, prefit$family$link)
    prefit_list[[key]] <- prefit
  }
  return(prefit_list)
}
