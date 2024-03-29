% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{readRDS_helper}
\alias{readRDS_helper}
\title{\code{base::readRDS} helper function that enables loading files from urls}
\usage{
readRDS_helper(file, .url_bool = FALSE, ...)
}
\arguments{
\item{file}{a \link[base]{connection} or the name of the file where the \R object
    is saved to or read from.}

\item{.url_bool}{bool whether the location is a url or not}

\item{...}{
  Arguments passed on to \code{\link[base:readRDS]{base::readRDS}}
  \describe{
    \item{\code{refhook}}{a hook function for handling reference objects.}
  }}
}
\description{
\code{base::readRDS} helper function that enables loading files from urls
}
