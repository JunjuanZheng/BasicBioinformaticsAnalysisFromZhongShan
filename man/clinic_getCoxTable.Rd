
\name{getCoxTable}
\alias{getCoxTable}
\title{get a model result of univariate and multivariate Cox regression Models}
\description{
 getCoxTable help get a model result of univariate and multivariate Cox regression Models.
}
\usage{
getCoxTable(data,
            time.col,
            status.col,
            stepwise=F,direction="both",
            cluster,control,
            dig=2,
            names = "test1")
}
\arguments{
  \item{data}{a data frame containing cluster,time and status value.}
  \item{time.col}{colnames of time value}
  \item{status.col}{colnames of status value}
  \item{stepwise}{Whether use stepwise strategy to simple multiple cox model}
  \item{direction}{method of stepwise.See MASS::stepAIC.}
  \item{cluster}{the cluster values like age and sex.}
  \item{dig}{2.the decimal place of output numeric value.}
  \item{names}{part of file name.}
}
\details{

}
\value{
a list containing results of univariate and multivariate models.
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

library(survival)
data(lung)
cluster = c("age","sex","ph.ecog","ph.karno","pat.karno")
## Not use stepwise
l1 <- getCoxTable(data=lung,
                  time.col="time",
                  status.col = "status",
                  cluster=cluster,
                  dig=2,
                  names = "a1")

## Use stepwise
l1 <- getCoxTable(data=lung,
                  time.col="time",
                  status.col = "status",
                  cluster=cluster,
                  stepwise=T,
                  dig=2,
                  names = "a2")
}
\keyword{getCoxTable}




