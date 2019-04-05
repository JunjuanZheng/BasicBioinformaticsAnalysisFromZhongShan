
\name{BalanceSplit}
\alias{BalanceSplit}
\title{split a data frame into balanced parts according to specified cluster factor}
\description{
 BalanceSplit help quickly splitting a data frame into balanced parts according to specified cluster factor.
}
\usage{
  BalanceSplit(data,
             cluster,
             nsplit = 2,
             seed.range=1:500,
             p.cutoff=0.1)
   summary.BalanceSplit(BalanceList,
                        upset.seed=c(3,11),
                        set=1)
}
\arguments{
  \item{data}{a data frame containing cluster cols.}
  \item{cluster}{cluster colnames.}
  \item{nsplit}{the number of balance splitting.Default is 2.}
  \item{seed.range}{candidate seeds.Default is 1:500.}
  \item{p.cutoff}{the cut off of P value.Default is 0.1.}
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
library(lucky)
Plus.library(c("tableone","knitr"))
data(data.TableOne)
data=data.TableOne
cluster = c("age","gender","his1","his2","his.grade","Stage","pT","N.status","M.status","chem.status")
L <- BalanceSplit(data,
                  cluster,
                  nsplit = 2,
                  seed.range=1:500,
                  p.cutoff=0.1)
table1 <- L[["seed_456"]][["table"]]
kable(table1[,1:3],align = 'c')
data.TableOne2 <- L[["seed_456"]][["data"]]

## use upset plot to visualize the difference between seeds
summary.BalanceSplit(BalanceList = L,
                     upset.seed=c(3,11),
                     set=1)
}
\keyword{BalanceSplit}




