% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moveData.R
\name{moveData}
\alias{moveData}
\title{Constructor of \code{moveData} objects}
\usage{
moveData(data)
}
\arguments{
\item{data}{A data frame with columns: \code{ID} (track ID(s)), \code{step}
(step length), \code{angle} (turning angle, if any), \code{x} (Easting or
longitude), \code{y} (Norting or latitude), and covariates (if any).}
}
\value{
An object of class \code{moveData}.
}
\description{
This constructor is not intended to be used, except inside the function
\code{\link{prepData}}. Refer to the documentation for that function.
}
