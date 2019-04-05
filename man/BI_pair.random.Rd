
\name{pair.random}
\alias{pair.random}
\title{pair.random}
\description{
 a simple method to make random-matched pair data by giving a contrast.
}
\usage{
  pair.random(data,
              contrast,
              cluster,
              seed = 2018,
              ratio=0.2,
              pair = TRUE,
              summary = T)
}
\arguments{
  \item{contrast}{string.colnames.}
  \item{seed}{seed used in sample()}
  \item{pair}{logic.If T,give a random-matched pair data by giving a contrast}
  \item{summary}{logic.If T,give a simple data.}

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
data(data.selfprocess2.2)
colnames(data.selfprocess2.2)
pr1 <- pair.random(data=data.selfprocess2.2,
                   contrast="N.status",
                   cluster=c("T.status","gender.status","age.status"),
                   seed = 2018,
                   ratio=0.5,
                   pair = TRUE,
                   summary = T)
View(pr1)
}
\keyword{pair.random}

