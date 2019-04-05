
\name{ensemblGoEnrich}
\alias{ensemblGoEnrich}
\title{Go annotation based on lucky::ensemblGo}
\description{
   ensemblGoEnrich help do a simple Go annotation for given genes based on ensemblGo
}
\usage{
  ensemblGoEnrich(genes,
                  parallel = F,
                  save.file = F,
                  names = "love")
}
\arguments{
  \item{genes}{a character/factor vector of genes}
  \item{parallel}{whether use parallel strategy to accelerate the process}
  \item{save.file}{whether save the data}
  \item{names}{part of saved file name}
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
Internet must be available for work
}
\seealso{
lucky::ensemblGo
}
\examples{
## test dataset
test <- common.annot$ENSEMBL;genes <- test[1:1000]

## Go enrich annotation
system.time(
  b <- ensemblGoEnrich(genes,parallel = F)
)

## View the result
View(b$merge)
View(b$metadata)
}

\keyword{ensemblGoEnrich}




