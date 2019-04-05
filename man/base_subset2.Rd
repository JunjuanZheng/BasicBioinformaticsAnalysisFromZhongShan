
\name{subset2}
\alias{subset2}
\title{get a subset of a data frame}
\description{
 subset2() can get a subset of a data frame,and keep data frame style regardless of the nrow or ncol of this subset.
}
\usage{
  subset2(data,
          colnames = NULL,
          rownames = NULL)
}
\arguments{
  \item{data}{a data frame}
  \item{colnames}{the colnames vector of a subset}
  \item{rownames}{the rownames vector of a subset}
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
## example
data <- data.frame(a=1:10,b=11:20,c=21:30,
                   row.names = paste("A",0:9,sep = ""))
subset2(data)

c = c("a","b");r=c("A0","A1","A2")
subset2(data,c,r)

c = c("a");r=c("A0","A1","A2")
subset2(data,c,r)

c = c("a","b");r=c("A0")
subset2(data,c,r)

c = c("a");r=c("A0")
subset2(data,c,r)
}
\keyword{subset2}




