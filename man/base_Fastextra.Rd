
\name{Fastextra}
\alias{Fastextra}
\title{a fast way to extra part of elements in a vetor.}
\description{
 a fast way to extra part of elements in a vetor.
}
\usage{
  Fastextra(vt,split,n)
}
\arguments{
  \item{vt}{a vector with common split}
  \item{split}{the common split}
  \item{n}{which splited string should be extrated}
}
\details{
a fast way to extra part of elements in a vetor.
}
\value{
 a vector with lots of extrated strings.
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
a <- rep(paste("a","b",sep="-"),10)
Fastextra(a,"-",1) # "a" "a" "a" "a" "a" "a" "a" "a" "a" "a"
Fastextra(a,"-",2) # "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"
}
\keyword{Fastextra}




