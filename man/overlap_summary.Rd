% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/range_overlap_join.R
\name{overlap_summary}
\alias{overlap_summary}
\title{Summarise events overlapping of 2 event data frames.}
\usage{
overlap_summary(x, y, by, sstable = FALSE, flextable = FALSE)
}
\arguments{
\item{x, y}{data.frames.}

\item{sstable}{logical value specifying whether to return in sstable format. Default is FALSE.}

\item{flextable}{logical value specifying whether to build flextable object. Default it FALSE. Set to TRUE forces sstable to TRUE.}

\item{ids}{the id col(s) wrapped in dplyr::vars() distinguishing different events.
Based on this is the summarisation taken.
If missing, except for the last 2, every columns will be automatically included into ids.}

\item{cols, cols.x, cols.y}{variables wrapped in dplyr::vars(), passed to overlap_join.}
}
\value{
A summary table that has at least 6 columns.
\itemize{
\item id column(s): identification of each event.
\item X: total time when X happens, regardless Y
\item Y: total time when Y happens, regardless X
\item X_Y: time duration when X and Y happen simultaneously.
\item X_notY: time duration when only X happens
\item notX_Y: time duration when only Y happens
}
}
\description{
A function to summarise the overlapping status of 2 event data frames
}
\seealso{
\link[dplyr]{vars}, \link[data.table]{foverlaps}, \link{overlap_join}
}
