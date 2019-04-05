
\name{Plus.order}
\alias{Plus.order}
\title{run order() for lots of cols in a matrix}
\description{
run order() for lots of cols in a matrix
}
\usage{
  Plus.order(data,col,decreasing=T)
}
\arguments{
  \item{data}{a matrix or data frame}
  \item{col}{colnames of the data}
}
\details{
run order() for lots of cols in a matrix
}
\value{
a matrix
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
data(fpkm)
col = c("patient25","patient8")
o = Plus.order(data,col,decreasing=T);View(o)
}
\keyword{Plus.order}
