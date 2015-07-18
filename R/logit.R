
#' Logit function
logit <- function(p)
{
  return(log(p/(1-p)))
}

#' Inverse logit function
inv.logit <- function(x)
{
  return(exp(x)/(1+exp(x)))
}
