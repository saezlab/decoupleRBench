% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run_benchmark}
\alias{run_benchmark}
\title{Benchmark pipeline built on the statistical method wrapper \link{decouple}.}
\usage{
run_benchmark(
  .design,
  .form = TRUE,
  .perform = TRUE,
  .minsize = 5,
  .silent = TRUE,
  .downsample_pr = FALSE,
  .downsample_roc = FALSE,
  .downsample_times = 100,
  .url_bool = FALSE
)
}
\arguments{
\item{.form}{bool whether to format or not}

\item{.perform}{bool whether to calculate ROC and performance summary}

\item{.minsize}{minimum size of each set}

\item{.silent}{bool whether to silence wanring messages}

\item{.downsample_pr}{whether to downsample precision recall curve TNs}

\item{.downsample_roc}{whether to downsample ROC true negatives}

\item{.downsample_times}{downsampling iterations}

\item{.url_bool}{bool whether the location is a url or not}
}
\value{
An S4 object of \link{BenchResult-class}
}
\description{
Benchmark pipeline built on the statistical method wrapper \link{decouple}.
}
\seealso{
See \link{input_tibble} for a description of the params/columns
of .design (i.e. input tibble).
}
