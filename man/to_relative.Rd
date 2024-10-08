% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_relative_time.R
\name{to_relative}
\alias{to_relative}
\alias{to_relative.Date}
\alias{to_relative.POSIXt}
\title{Convert from absolute date/time to relative numeric number}
\usage{
to_relative(..., zero = min(...), shift = 0, units = "auto")

\method{to_relative}{Date}(..., zero = min(...), shift = 0)

\method{to_relative}{POSIXt}(..., zero = min(...), shift = 0)
}
\arguments{
\item{...}{variables of type Date/POSIXt}

\item{zero}{A value specifying the date or time zero}

\item{shift}{If zero is not date/time 0, but is forward/backward-shifted. Changing the shift function can correct that.}

\item{units}{A character vector passed to \link{difftime}, default is auto}
}
\value{
Vectors of type double. Number of outputted vectors is equal to number of inputted ones
}
\description{
A function to convert calendar date/time to relative day/time(s). Useful in survival analyses
}
