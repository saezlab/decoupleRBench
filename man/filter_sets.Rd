% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{filter_sets}
\alias{filter_sets}
\title{Helper Function to filter and format the gene set resource}
\usage{
filter_sets(set_source, source_col, filter_col, filter_crit, .minsize, .silent)
}
\arguments{
\item{set_source}{Set Source (e.g. TF regulon sets, GO:term sets, etc)}

\item{source_col}{name of the column with the source for the set source}

\item{filter_col}{name of the column by which we wish to filter}

\item{filter_crit}{criteria by which we wish to filter the \code{filter_col}}

\item{.minsize}{minimum size of each set}

\item{.silent}{bool whether to silence wanring messages}
}
\value{
returns a filtered and formatted set source
}
\description{
Helper Function to filter and format the gene set resource
}
\details{
Filtering can be omitted if \code{filter_col} is \code{NA}.
}
