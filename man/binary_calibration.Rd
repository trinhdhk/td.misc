% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{binary_calibration}
\alias{binary_calibration}
\title{Calibration curve for binary outcome}
\usage{
binary_calibration(
  pred,
  pred_rep = NULL,
  obs,
  title = NULL,
  method = c("loess", "splines"),
  se = TRUE,
  hist = TRUE,
  hist.normalize = TRUE,
  yscale = 0.1,
  ...,
  theme = ggplot2::theme_bw
)
}
\arguments{
\item{pred}{predictor values/probabilities}

\item{pred_rep}{for K-Fold CV, plot the predicted values for all repetitions}

\item{obs}{observed values}

\item{method}{either loess or splines}

\item{se}{whether to draw se ribbon}

\item{hist}{whether to draw histograms for each class}

\item{hist.normalize}{whether to normalize the height of histograms to make them reflects the true unbalance of the two class.}

\item{...}{control parameters passed to the smoothing method, commonly span and knots}

\item{theme}{theme function. default to theme_bw}
}
\value{
a ggplot
}
\description{
Helper Function to create calibration curve for binary outcomes
}
