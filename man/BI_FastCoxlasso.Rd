
\name{FastCoxlasso}
\alias{FastCoxlasso}
\title{Penalized Cox regression model via lasso algrithm}
\description{
    FastCoxlasso help do a plot-based penalized Cox regression model via lasso algrithm.
}
\usage{
FastCoxlasso(
  expr.matrix,
  design,
  select,
  event.status=c("TTR.status","DFS.status","OS.status")[2],
  event.time=c("TTR.time","DFS.time","OS.time")[2],
  event.lower =c(89,89,89),
  lambda.n = 3,ln.lambda = 0,label=F,
  digits = 5,
  names = "love"
)
}
\arguments{
\item{expr.matrix}{expression matrix}
\item{design}{design object}
\item{select}{selected markers or genes}
\item{event.status}{colnames of event status}
\item{event.time}{colnames of event time}
\item{event.lower}{the lower limit of event time}
\item{lambda.n}{the counts of genes with parameters you want}
\item{ln.lambda}{ln(lambda)}
\item{label}{whether to show the genes label on the plot.Default is F}
\item{digits}{the digits of model coefficient}
\item{names}{part of saved files name}
}
\details{

}
\value{
 a list containing gene-based survival data and model data
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
expr.matrix = train.fpkm
design = design.train
select = rocgenes
event.status="DFS.status"
event.time="DFS.time"
event.lower = 89
label = T
label=TRUE
digits = 5
names = "love"
CO1 <- FastCoxlasso(expr.matrix,
                    design,
                    select,
                    event.status,
                    event.time,
                    event.lower,
                    lambda.n = 5,
                    ln.lambda= -3.55,
                    label,digits,names)
View(CO1$DFS$modeldata)

## Common Use
FC1 <- FastCoxlasso(expr.matrix = train.fpkm,
                    design= design.train,
                    select = rocgenes,
                    event.status="DFS.status",
                    event.time="DFS.time",
                    event.lower= 89,
                    lambda.n = 5,
                    ln.lambda= -3.55,
                    names = "train")

}
\keyword{FastCoxlasso}




