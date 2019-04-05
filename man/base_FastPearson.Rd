
\name{FastPearson}
\alias{FastPearson}
\title{Get correlation index from expression matrix}
\description{
FastPearson help get cor from expression matrix.
}
\usage{
FastPearson(data,
            transposition = F,
            control.markers=NULL,
            target.markers=NULL,
            method = "pearson",
            order = F,
            parallel = F)
}
\arguments{
\item{data}{a expression matrix}
\item{transposition}{Whether transpose the matrix.IF you input a gene expression matrix with rownames gene,you should set transposition = T}
\item{control.markers}{If NULL,use all pairs strategy,else use one pair strategy.Only one control marker is supported}
\item{target.markers}{when control.marker is not a NULL value,it specified the target variates you want.Default is NULL,which means you wan correlation between control.markers to other variates}
\item{method}{One of "pearson","kendall" and "spearman".Default is "pearson"}
\item{order}{whether order the cormatrix with corrplot::corrMatOrder}
\item{parallel}{whether use parallel strategy}
}
\details{

}
\value{
 a list of long data frame and a cormatrix
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

library(lucky)
data=state.x77;colnames(data)

## All pairs
l1 <- FastPearson(data = state.x77,
            transposition = F,
            control.markers=NULL,
            target.markers=NULL,
            method = "pearson",
            order = F,
            parallel = F)
View(l1$longdata)
View(l1$cormatrix)

## One pair with no specified target
l2 <- FastPearson(data = state.x77,
                  transposition = F,
                  control.markers="Income",
                  target.markers=NULL,
                  method = "pearson",
                  order = F,
                  parallel = F)
View(l2$longdata)

## One pair with specified target
l3 <- FastPearson(data = state.x77,
                  transposition = F,
                  control.markers="Income",
                  target.markers=c("Illiteracy","Life Exp"),
                  method = "pearson",
                  order = F,
                  parallel = F)
View(l3$longdata)

}
\keyword{FastPearson}




