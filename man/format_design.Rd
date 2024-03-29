% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{format_design}
\alias{format_design}
\title{Helper Function to to generate the bools used to check if the current
locations/rds objects are the same as the previous one.}
\usage{
format_design(.design)
}
\arguments{
\item{.design}{input tibble used to provide the experimental design for each
benchmark run}
}
\description{
Helper Function to to generate the bools used to check if the current
locations/rds objects are the same as the previous one.
}
\details{
This is used to limit the number of times that any of the
prerequsites is loaded.
}
\keyword{internal}
