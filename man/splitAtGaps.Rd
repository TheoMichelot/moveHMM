% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitAtGaps.R
\name{splitAtGaps}
\alias{splitAtGaps}
\title{Split track at gaps}
\usage{
splitAtGaps(data, maxGap = 60, shortestTrack = 0, units = "mins")
}
\arguments{
\item{data}{Data frame with (at least) columns for "ID" and "time"}

\item{maxGap}{Longest allowed gap, in minutes by default (but see argument
"units"). Track will be split at longer gaps.}

\item{shortestTrack}{Shortest track to keep after splitting, in minutes.
Shorter tracks will be removed from the output data set.}

\item{units}{Character string, e.g., "mins" (default), "secs", "hours",
"days". Time units used for the maxGap and shortestTrack arguments.}
}
\value{
Data frame with identical structure as input, where ID column
has been replaced by new ID for split tracks. Old ID still accessible as
ID_old column
}
\description{
Defines a new ID variable, which changes where there is a long gap in
the data. This is sometimes useful for preprocessing data prior to using
\code{prepData}.
}
