
\name{getOneDesignEset}
\alias{getOneDesignEset}
\title{get one DesignEset}
\description{
 get one DesignEset when one or lots of conditions exist in DesignEset.
}
\usage{
  getOneDesignEset(DesignEset,contrast.i)
}
\arguments{
\item{DesignEset}{a result from AddDesignEset()}
\item{n}{which condition and related DesignEset you want to get}
}
\details{
 getOneDesignEset help to get a DesignEset with only one condition.
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
 AddDesignEset()
}
\examples{
 DesignEset.i <- getOneDesignEset(DesignEset,c("N1","N0"))
}
\keyword{getOneDesignEset}

