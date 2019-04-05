
\name{FPKMtoTPM}
\alias{FPKMtoTPM}
\title{FPKM value converted to TPM value}
\description{
  FPKM value converted to TPM value
}
\usage{
  FPKMtoTPM(mt)
}
\arguments{
  \item{mt}{data frame or matrix.It always be a gene expression matrix.}
}
\details{

}
\value{
 a matrix contain TPM Values
}
\references{

}
\author{
  Weibin Huang
}
\note{

}
\seealso{
 data(fpkm)
}
\examples{
data(fpkm) # fpkm <- assay(data)
tpm <- FPKMtoTPM(fpkm)
tpm[1,1] == (fpkm[1,1]/sum(as.numeric(as.character(fpkm[,1]))))*(10^6) ## T
colSums(tpm)/nrow(tpm) # a value that should be a repeats of constant values.

}
\keyword{FPKMtoTPM}

