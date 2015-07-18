
#' Logit function
#' Used in n2w.
logit <- function(p)
{
  return(log(p/(1-p)))
}

#' Inverse logit function
#' Used in w2n.
inv.logit <- function(x)
{
  return(exp(x)/(1+exp(x)))
}
