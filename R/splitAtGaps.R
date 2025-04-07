
#' Split track at gaps
#'
#' Defines a new ID variable, which changes where there is a long gap in
#' the data. This is sometimes useful for preprocessing data prior to using
#' \code{prepData}.
#'
#' @param data Data frame with (at least) columns for "ID" and "time"
#' @param maxGap Longest allowed gap, in minutes by default (but see argument
#' "units"). Track will be split at longer gaps.
#' @param shortestTrack Shortest track to keep after splitting, in minutes.
#' Shorter tracks will be removed from the output data set.
#' @param units Character string, e.g., "mins" (default), "secs", "hours",
#' "days". Time units used for the maxGap and shortestTrack arguments.
#'
#' @return Data frame with identical structure as input, where ID column
#' has been replaced by new ID for split tracks. Old ID still accessible as
#' ID_old column
#'
#' @export
splitAtGaps <- function(data, maxGap = 60, shortestTrack = 0, units = "mins") {
    # Number of tracks
    nTracks <- length(unique(data$ID))

    # Save old ID and reinitialise ID column
    data$ID_old <- data$ID
    data$ID <- character(nrow(data))

    # Loop over tracks (i.e., over IDs)
    for(i_track in 1:nTracks) {
        # Indices for this track
        indThisTrack <- which(data$ID_old == unique(data$ID_old)[i_track])
        trackLength <- length(indThisTrack)

        # Time intervals
        dtimes <- difftime(data$time[indThisTrack[-1]],
                           data$time[indThisTrack[-trackLength]],
                           units = units)

        # Indices of gaps longer than maxGap
        indGap <- c(0, which(dtimes > maxGap), trackLength)

        # Create new ID based on split track
        subtrack_ID <- rep(1:(length(indGap) - 1), diff(indGap))
        data$ID[indThisTrack] <- paste0(data$ID_old[indThisTrack],
                                        "-", subtrack_ID)
    }

    # Only keep sub-tracks longer than some duration
    trackLengths <- sapply(unique(data$ID), function(id) {
        ind <- which(data$ID == id)
        difftime(data$time[ind[length(ind)]],
                 data$time[ind[1]],
                 units = units)
    })
    ID_keep <- names(trackLengths)[which(trackLengths >= shortestTrack)]
    data <- data[which(data$ID %in% ID_keep),]

    return(data)
}
