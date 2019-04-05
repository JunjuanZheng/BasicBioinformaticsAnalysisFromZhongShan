
\name{clinic.filters}
\alias{clinic}
\title{clinic.filters}
\description{
 Filtering clinical data.
}
\usage{
  clinic.filters(data,filter,
                filter.parameter,
                summary = T)
}
\arguments{
  \item{data}{a data frame}
  \item{filter}{string.names of columns of a numeric or character value in the data.}
  \item{filter.parameter}{strings containing operational characters such as ">","<","="}
  \item{summary}{logic.If T,print simple data after filtering,else complete data.}
}
\details{

}
\value{

}
\references{

}
\author{
  Weibin Huang
}
\note{

}
\seealso{

}
\examples{
 data(data.selfprocess2)
colnames(data.selfprocess2)
data.selfprocess2$LN.sum <- as.numeric(as.character(data.selfprocess2$LN.sum))
data.selfprocess2$OS.time <- as.numeric(as.character(data.selfprocess2$OS.time))
cf1 <- clinic.filters(data = data.selfprocess2,
                      filter = c("LN.sum","OS.time"),
                      filter.parameter = c(">15",">179"),
                      summary = F)
cf2 <- clinic.filters(data = data.selfprocess2,
                      filter = c("LN.sum","OS.time","BMI"),
                      filter.parameter = c(">15",">179",">15"),
                      summary = F)
View(cf1);View(cf2)
}
\keyword{clinic.filters}

