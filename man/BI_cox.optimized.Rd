
\name{cox.optimized}
\alias{cox.optimized}
\title{help get a best model based on prognosis information via Cox model}
\description{
help get a best model based on prognosis information via Cox model
}
\usage{
cox.optimized(
  expr.matrix,
  design,
  select,
  event.status=c("TTR.status","DFS.status","OS.status"),
  event.time=c("TTR.time","DFS.time","OS.time"),
  event.lower =c(89,89,89),
  direction="both",
  method = "pca",
  digits = 5
)
}
\arguments{
\item{expr.matrix}{expression matrix}
\item{design}{design object}
\item{select}{selected markers or genes}
\item{event.status}{colnames of event status}
\item{event.time}{colnames of event time}
\item{event.lower}{the lower limit of event time}
\item{direction}{parameter of MASS::stepAIC.One of "both", "backward", "forward".Default is "both"}
\item{method}{whether use PCA strategy to avoid multicollinearity.Default is "pca"}
\item{digits}{the digits of Cor}
}
\details{

}
\value{

}
\references{

}
\author{
  Weibin Huang<654751191@qq.com>
}
\note{

}
\seealso{

}
\examples{

}
\keyword{cox.optimized}




