
\name{cox.optimized4}
\alias{cox.optimized4}
\title{get survival significant markers based on KM-Curve data}
\description{
cox.optimized4 help get survival significant markers based on KM-Curve data
}
\usage{
cox.optimized4(expr.matrix,
               design,
               time.col,
               status.col,
               strategy  = c("median","maxstat")[1],
               smethod="LogRank",
               pmethod="HL",
               cutoff.p = 0.05,
               digits = 5,
               parallel = F)
}
\arguments{
\item{expr.matrix}{expression matrix}
\item{design}{design object}
\item{time.col}{colname of time value}
\item{status.col}{colname of status value}
\item{strategy}{one of "median" and "maxstat"}
\item{smethod}{Default is "LogRank"}
\item{pmethod}{Default is "HL"}
\item{cutoff.p}{p value cut-off}
\item{digits}{digits of result}
\item{parallel}{whether to use parallel strategy when strategy is median.Note that the parallel would always be applied with maxstat method}
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
the maxstat function depends on maxstat package.Parameters like smethod and pmethod are belongs to maxstat function.
}
\seealso{
maxstat::maxstat.test

}
\examples{
## This is a simulative process and NOT RUN
cox.optimized4(expr.matrix,
               design,
               time.col,
               status.col,
               strategy  = "median",
               cutoff.p = 0.05,
               digits = 5,
               parallel = T)

## This is a simulative process and NOT RUN
cox.optimized4(expr.matrix,
               design,
               time.col,
               status.col,
               strategy  = "maxstat",
               smethod="LogRank",
               pmethod="HL",
               cutoff.p = 0.05,
               digits = 5)

}
\keyword{cox.optimized4}

