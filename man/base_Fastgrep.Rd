
\name{Fastgrep}
\alias{Fastgrep}
\title{A Fast grep when lots of pattern used}
\description{
 A Fast grep when lots of pattern used
}
\usage{
  Fastmatch(pattern,x)
  Fastgrep(pattern,x)
}
\arguments{
  \item{pattern}{a string}
  \item{x}{a vector that may contain patterns.}
}
\details{
A Fast grep when lots of pattern used.return an unique position of all the pattern.Fastgrep() is based on base::grep().
}
\value{
a numeric of position
}
\references{
grep()
}
\author{
  Weibin Huang
}
\note{

}
\seealso{
grep()
}
\examples{
## compared with grep and Fastgrep
grep("a",c("b","c","d"))
Fastgrep("a",c("b","c","d"))
grep("b",c("b","c","d"))
Fastgrep("b",c("b","c","d"))
Fastgrep(c("a","c"),c("b","c","d"))
Fastgrep(c("b","c"),c("b","c","d"))
}
\keyword{Fastgrep()}




