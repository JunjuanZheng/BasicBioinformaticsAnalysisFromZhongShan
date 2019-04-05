
\name{cox.optimized3}
\alias{cox.optimized3}
\title{select prognosis-predictive genes}
\description{
    cox.optimized3 help select prognosis-predictive genes base on expression matrix and related design object
}
\usage{
cox.optimized3(expr.matrix,
               design,
               sd.cutoff = NULL,
               p.cutoff = 0.0001,
               FDR.cutoff = NULL,
               event.status=c("TTR.status","DFS.status","OS.status"),
               event.time=c("TTR.time","DFS.time","OS.time"),
               event.lower =c(89,89,89),
               show.music =T)
}
\arguments{
  \item{expr.matrix}{expression matrix}
  \item{design}{design object}
  sd.cutoff
  \item{sd.cutoff}{the cut-off of standard deviation.Default is NULL,and you can set a number if you want to use a sd filter}
  \item{p.cutoff}{the cut-off of P value.Default is 0.0001}
  \item{FDR.cutoff}{the cut-off of BH adjusted P value.Default is 0.001.If NULL,the FDR filter would not be applied}
  \item{event.status}{character, the colnames of event status}
  \item{event.time}{]character.the colnames of event time}
  \item{event.lower}{numeric.the lower border of time}
  \item{show.music}{whether show music at the end of the programe}
}
\details{
cox.optimized3 depends on survival package and some parrallel algrithm
}
\value{
a list containing cox model genes and related cox result
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
## This is a simulative process and NOT RUN
cox.model <- cox.optimized3(expr.matrix,
                            design,
                            p.cutoff = 0.0001,
                            FDR.cutoff = 0.001,
                            event.status=c("TTR.status","DFS.status","OS.status"),
                            event.time=c("TTR.time","DFS.time","OS.time"),
                            event.lower =c(89,89,89))
save(cox.model,file = "seed35_cox.model.rda")

## Garbage management
garbage = c("expr.matrix1","design1","event.status.i","event.time.i","event.lower.i")
for(n in garbage){rm(list = ls(pattern = n),envir = .GlobalEnv)}
}
\keyword{cox.optimized3}




