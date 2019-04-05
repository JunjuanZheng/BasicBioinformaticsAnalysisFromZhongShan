
\name{stringposition}
\alias{stringposition}
\title{find the position of a specified string pattern}
\description{
stringposition help find the position of a specified string pattern
}
\usage{
  stringposition(pattern,
               x,
               unlist=T)
}
\arguments{

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
If there is no specified pattern in a string,NA would be returned.
}
\seealso{

}
\examples{
stringposition(c("a","-"),
               "abc-bace-a")#unlist=T
stringposition(c("a","-"),
               "abc-bace-a",
               unlist = F)#unlist=F
stringposition("-", "abc")#
stringposition(c("d","-"), "abc-bace-a")#-中的位置
}
\keyword{stringposition}




