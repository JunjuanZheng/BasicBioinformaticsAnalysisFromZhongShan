
\name{QC.cluster}
\alias{QC.cluster}
\title{QC.cluster}
\description{
 QC.cluster()通过选取少数基因进行聚类，从而判断分组的差异程度。QC.cluster()基于hwb_BI::heatmap.dds()函数。输出结果为heatmap图。
}
\usage{
QC.cluster(DesignEset,
           contrast,
           select.c=1000,#选取的基因数
           seed=2018,
           savefile=F,
           log.convert = F,
           names="GSE1231",#part of file name
           cluster_rows = T,
           cluster_cols =T,
           row.k = 4)
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
QC.cluster(DesignEset,
           contrast,
           select.c=1000,#选取的基因数
           seed=2018,
           savefile=F,
           log.convert = F,
           names="GSE1231",#part of file name
           cluster_rows = T,
           cluster_cols =F,
           row.k = 4)
}
\keyword{QC.cluster}

