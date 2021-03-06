% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_classify.R
\name{classify}
\alias{classify}
\title{Classify a vector by specified classifier}
\usage{
classify(vector, classifier)
}
\arguments{
\item{vector}{a vector}

\item{classifier}{a list of classifier}
}
\description{
classify help classify a vector by specified classifier
}
\examples{
## example
vector <- 1:100

## classifier exist
classifier1 <- list(
  lower = c(1:49),
  upper = c(50:90)
)
v1 <- classify(vector,classifier1)
table(v1)

## upper not exist
classifier2 <- list(
  lower = c(1:49),
  upper = c(101:900)
)
v1 <- classify(vector,classifier2)
table(v1)
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
