
\name{cox.threshold3}
\alias{cox.threshold3}
\title{Get a comprimized cut-off for survival analysis among mutiple design object}
\description{
Get a comprimized cut-off for survival analysis among mutiple design object
}
\usage{
cox.threshold3(expr.matrix=data.fpkm,
               design.list,
               model,
               event.status=status,
               event.time=time,
               event.lower = event.lower,
               interval.n  = 1000,
               cutoff.pval = 0.05)
}
\arguments{
\item{expr.matrix}{expression object}
\item{design.list}{a list of design object.The names of them must be provided}
\item{model}{model from List.ModelFromOptimized2}
\item{event.status}{colname of status}
\item{event.time}{colname of time}
\item{event.lower}{lower limit of time}
\item{interval.n}{the cut-off selected parameter.Default is 1000}
\item{cutoff.pval}{cut-off of p value}
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

## These are simulative process and NOT RUN

## Quick Start
ct3 <- cox.threshold3(expr.matrix=data.fpkm,
                      design.list = list(train = design.train,
                                         valid = design.valid),
                      model = COX2.2$DFS,
                      event.status=status,
                      event.time=time,
                      event.lower = event.lower,
                      interval.n  = 1000,
                      cutoff.pval = 0.05)

## parameters
model =COX2.2$DFS
expr.matrix=data.fpkm
design.list = list(train = design.train,valid = design.valid)
event.status=status
event.time=time
event.lower = event.lower
interval.n  = 1000
cutoff.pval = 0.05

}
\keyword{}




