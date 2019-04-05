
\name{QuickStart}
\alias{QuickStart}
\title{Run series R script automatically}
\description{
QuickStart help run series R script automatically
}
\usage{
QuickStart(script.path,
           data.path,
           result.path)
}
\arguments{
  \item{script.path}{the path of script}
  \item{data.path}{the path of data}
  \item{result.path}{the path of result}
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

## Delete median variables
list = ls()
keep = c("script.path","data.path","result.path.all","List.MODEL","names","seeds")
list1 = setdiff(list,keep)
rm(list = list1)
}
\keyword{QuickStart}




