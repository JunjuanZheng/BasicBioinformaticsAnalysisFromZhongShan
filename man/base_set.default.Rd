
\name{set.default}
\alias{set.default}
\title{help self-defined function to set default values}
\description{
 help self-defined function to set default values
}
\usage{
  set.default(input,default,output)
}
\arguments{
 \item{input}{a list contain input value}
 \item{default}{a list contain defalut value}
 \item{output}{a list contain value that you want to output when the input was as the same as default.}
}
\details{

}
\value{

}
\references{

}
\author{
  Weibin Huang
}
\note{

}
\seealso{
  Fastcorrplot()
}
\examples{
# example :set default for Fastcorrplot()
default <- c(rep(list(NULL),7))
library(grDevices);col <- colorRampPalette(c("blue","white","red"))(100)
input <- list(lower.col = lower.col,upper.col =upper.col,upper = upper,tl.pos = tl.pos,tl.col=tl.col,tl.srt=tl.srt,diag = diag)
do <- list(lower.col = col,upper.col =col,upper = "ellipse",tl.pos = "lt",tl.col="black",tl.srt=45,diag = "l")
output <- set.default(input,default,output)
}
\keyword{is.all.na}




