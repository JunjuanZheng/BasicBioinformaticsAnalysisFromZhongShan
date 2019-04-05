
\name{extra.surv}
\alias{extra.surv}
\title{Get survival data of a gene based on expression and design objects}
\description{
    extra.surv help get survival data of a gene based on expression and design objects.Note that multiple genes would not be supported.
}
\usage{
extra.surv(expr.matrix,
                 design,
                 select,
                 event.status,
                 event.time,
                 event.lower,
                 mode=c("median","maxstat")[1])
}
\arguments{
  \item{expr.matrix}{expression matrix}
  \item{design}{design object}
  \item{select}{a gene.Multiple genes are not supported}
  \item{event.status}{a colname of status}
  \item{event.time}{a colname of time}
  \item{event.lower}{the lower limit of the time}
  \item{mode}{Default is "median".You can also set a number of percentum or a decimal less than 1.When mode = "maxstat",the maxstat automatical selected cut-off strategy would be used}
}
\details{

}
\value{
    a list containing P value of significant in univariate cox model and the metadata of survival information and binary gene expression
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
expr.matrix=train.fpkm
design = design.train
select = rocgenes[1]
event.status = "DFS.status"
event.time = "DFS.time"
event.lower = 89
mode="median"

## method = "median"
l1 <- extra.surv(expr.matrix,
                 design,
                 select,
                 event.status,
                 event.time,
                 event.lower,
                 mode="median")

## method = "maxstat"
l1 <- extra.surv(expr.matrix,
                 design,
                 select,
                 event.status,
                 event.time,
                 event.lower,
                 mode="maxstat")
}
\keyword{extra.surv}




