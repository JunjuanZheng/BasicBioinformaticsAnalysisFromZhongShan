% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_convertStatus.R
\name{convertStatus}
\alias{convertStatus}
\title{convertStatus}
\usage{
convertStatus(data, time, status, event = "1", cutoff = c(365, 365 *
  2), names = c("Year_1", "Year_2"), filter = T)
}
\arguments{
\item{data}{the data containing time and status columns}

\item{time}{colname of time}

\item{status}{colname of status}

\item{event}{the symbol of event happening}

\item{cutoff}{the cut off of time to create a binary variable}

\item{filter}{if \code{filter = T},then the sample that with less time and event not happening would be consider as \code{NA},which means that we could not do a correct judge.}

\item{name}{the colnames of new binary variables}
}
\value{
a data frame of binary variables
}
\description{
convert time and status to a binary variable
}
\examples{
library(lucky)
data("rna.design")
df <- convertStatus(data = rna.design,
                    time = "OS.time",status = "OS.status",
                    event = "1",
                    cutoff = c(365,365*2),
                    names = c("Year_1","Year_2"))
View(df)
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
