
\name{result.ROCList}
\alias{result.ROCList}
\title{show the progonis value of multiple model genes from ROCList object of bootstrap.lasso function}
\description{
show the progonis value of multiple model genes from ROCList object of bootstrap.lasso function
}
\usage{
result.ROCList(expr.matrix.list,
               design.list,
               ROCList,
               contrast="N.status",
               contrast.control="N0",
               palette=c(1,4),
               save.file = T,
               width = 12,height=10,
               names = "love")
}
\arguments{
  \item{expr.matrix.list}{the list of comparing expression matrix}
  \item{design.list}{the list of design object of comparing expression matrix}
  \item{L1}{the result of lucky::bootstrap.lasso}
  \item{contrast}{the colnames of contrast}
  \item{contrast.control}{the names of control object in the contrast}
  \item{return.freq}{If F,return roc related plot;If T,return data}
  \item{palette}{a integer vector to select color from lucky::mycolor}
  \item{save.file}{whether save a PDF plot}
  \item{width,height}{the width and height of saved plot}
  \item{names}{part of file name}
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

}
\keyword{result.ROCList}




