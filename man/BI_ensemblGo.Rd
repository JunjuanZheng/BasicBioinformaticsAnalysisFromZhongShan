
\name{ensemblGo}
\alias{ensemblGo}
\title{Annotate the ENSG id with GO terms via ensembl online}
\description{
ensemblGo annotate the ENSEMBL ID with GO terms via ensembl online
}
\usage{
  ensemblGo(genes,
            parallel = F,
            save.file = F,
            names = "love")
}
\arguments{
  \item{genes}{a character vector of ENSG id}
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

}
\seealso{

}
\examples{

## test genes
test <- common.annot$ENSEMBL;genes <- test[1:50]

## Not use parallel
system.time(
  a <- ensemblGo(genes)
)
# User  System  time (50 tests)
# 0.56  0.15    36.20


## Use parallel
system.time(
  a <- ensemblGo(genes,parallel = T)
)
# User  System  time (50 tests)
# 0.09  0.01    14.19

## Get available Go annotation
a2 <- a[a$Accession != "NotAvailable",]

}
\keyword{ensemblGo}




