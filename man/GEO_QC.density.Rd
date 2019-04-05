
\name{QC.density}
\alias{QC.density}
\title{QC.density}
\description{
 QC.density可以观察不同分组（treat vs. control）的芯片表达量密度。基于ggpubr的ggdensity函数。
}
\usage{
QC.density(DesignEset,
           contrast,
           savefile=F,
           names="test1",
           limit.array = 30,
           limit.genes=3000,
           seed=2018)
}
\arguments{

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

}
\examples{
data(DesignEset)
contrast=list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")) # a list that contain lots of contrasts.
QC.density(DesignEset,
           contrast,
           savefile=F,
           names="test1",
           limit.array = 30,
           limit.genes=3000,
           seed=2018)
}
\keyword{QC.density}

