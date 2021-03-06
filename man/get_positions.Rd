% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liftover.r
\name{get_positions}
\alias{get_positions}
\title{Lookup positions for given rsids in particular build}
\usage{
get_positions(
  rsid,
  build = 37,
  method = c("opengwas", "biomart")[1],
  splitsize = 50000
)
}
\arguments{
\item{rsid}{rsid}

\item{build}{build (36, 37 default or 38)}

\item{method}{"opengwas" (fastest) or "biomart"}

\item{splitsize}{Default 50000}
}
\value{
data frame
}
\description{
Lookup positions for given rsids in particular build
}
