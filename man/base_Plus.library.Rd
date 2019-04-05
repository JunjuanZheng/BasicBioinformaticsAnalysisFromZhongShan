
\name{Plus.library}
\alias{Plus.library}
\title{Plus.library}
\description{
 A fast way for Chinese people to automatically install and library package.
}
\usage{
  Plus.library(packages)
}
\arguments{
  \item{packages}{a vector with package names}
}
\details{
  Plus.library() check localed installed packages.If you haven't install the package,it would install it automatically.If packages have installed,it libraries the packages.Plus.library depends on pacman packages.
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

}
\examples{
  Plus.library(c("DT","RColorBrewer","stringr"))
}
\keyword{Plus.library}

