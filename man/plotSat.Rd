% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSat.R
\name{plotSat}
\alias{plotSat}
\title{Plot observations on satellite image}
\usage{
plotSat(
  data,
  zoom = NULL,
  location = NULL,
  segments = TRUE,
  compact = TRUE,
  col = NULL,
  alpha = 1,
  size = 1,
  states = NULL,
  animals = NULL,
  ask = TRUE,
  return = FALSE
)
}
\arguments{
\item{data}{Data frame of the data, with necessary fields 'x' (longitude values) and
'y' (latitude values).}

\item{zoom}{The zoom level, as defined for \code{get_map}. Integer value between
3 (continent) and 21 (building).}

\item{location}{Location of the center of the map to be plotted.}

\item{segments}{\code{TRUE} if segments should be plotted between the observations (default),
\code{FALSE} otherwise.}

\item{compact}{\code{FALSE} if tracks should be plotted separately, \code{TRUE}
otherwise (default).}

\item{col}{Palette of colours to use for the dots and segments. If not specified, uses default palette.}

\item{alpha}{Transparency argument for \code{geom_point}.}

\item{size}{Size argument for \code{geom_point}.}

\item{states}{A sequence of integers, corresponding to the decoded states for these data
(such that the observations are colored by states).}

\item{animals}{Vector of indices or IDs of animals/tracks to be plotted.
Default: \code{NULL}; all animals are plotted.}

\item{ask}{If \code{TRUE}, the execution pauses between each plot.}

\item{return}{If \code{TRUE}, the function returns a ggplot object (which can be edited and
plotted manually). If \code{FALSE}, the function automatically plots the map (default).}
}
\description{
Plot tracking data on a satellite map. This function only works with longitude
and latitude values (not with UTM coordinates), and uses the package \code{ggmap}
to fetch a satellite image from Google. An Internet connection is required to use
this function.
}
\details{
If the plot displays the message "Sorry, we have no imagery here", try a
lower level of zoom.
}
\references{
D. Kahle and H. Wickham. ggmap: Spatial Visualization with ggplot2.
The R Journal, 5(1), 144-161.
URL: http://journal.r-project.org/archive/2013-1/kahle-wickham.pdf
}
