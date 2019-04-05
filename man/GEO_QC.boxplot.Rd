
\name{QC.boxplot}
\alias{QC.boxplot}
\title{QC.boxplot}
\description{
 质量控制相关的箱式图。limit.array = 30,limit.genes=3000指的是如果array个数>30，就不进行全基因分析，只随机选取3000个基因进行质量分析。
}
\usage{
QC.boxplot(DesignEset,
           contrast,
           savefile=F,
           names="test1",
           limit.array = 30,
           limit.genes=3000,
           seed=2018)
}
\arguments{
 \item{contrast}{list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")).a list that contain lots of contrasts}
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
QC.boxplot(DesignEset,
           contrast,
           savefile=F,
           names="test1",
           limit.array = 30,
           limit.genes=3000,
           seed=2018)
}
\keyword{QC.boxplot}

