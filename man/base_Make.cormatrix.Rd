
\name{Make.cormatrix}
\alias{Make.cormatrix}
\title{translation between cormatirx and related long data frame}
\description{
 translation between cormatirx and related long data frame
}
\usage{
  Make.cormatrix(data,
               source.col=NULL,
               target.col=NULL,
               value=NULL)
}
\arguments{
  \item{data}{a cormatrix or a long data frame from a cormatirx}
  \item{source.col}{If data is a cormatrix,it isn't needed.If data is a long data frame,it is the colnames of source columns.}
  \item{target.col}{If data is a cormatrix,it isn't needed.If data is a long data frame,it is the colnames of target columns.}
  \item{value}{If data is a cormatrix,it isn't needed.If data is a long data frame,it is the colnames of value columns.}
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
### example
## get cormatrix
data=state.x77;colnames(data)
result1 <- FastMINE(data,
                    transposition = F,
                    method = "all.pairs",
                    control.markers="Income",
                    target.markers=NULL)
data <- result1[["MIC.matirx"]]

## cormatrix to long data frame
data1 <- Make.cormatrix(data)
colnames(data1) #[1] "source" "target" "value"

## long data frame to cormatrix
data2 <- Make.cormatrix(data=data1,
                        source.col="source",
                        target.col="target",
                        value="value")
}
\keyword{Make.cormatrix}




