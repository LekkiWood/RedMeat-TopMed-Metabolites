get_context_effect <- function(model, pm_name, cwc_name) {
  
  b <- fixef(model)
  V <- vcov(model)
  
  beta_pm  <- b[pm_name]
  beta_cwc <- b[cwc_name]
  
  beta_context <- beta_pm - beta_cwc
  
  var_pm  <- V[pm_name,  pm_name]
  var_cwc <- V[cwc_name, cwc_name]
  covar   <- V[pm_name,  cwc_name]
  
  se_context <- sqrt(var_pm + var_cwc - 2 * covar)
  
  z  <- beta_context / se_context
  p  <- 2 * (1 - pnorm(abs(z)))
  
  ci_lower <- beta_context - 1.96 * se_context
  ci_upper <- beta_context + 1.96 * se_context
  
  data.frame(
    beta_context = beta_context,
    se = se_context,
    z = z,
    p = p,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}