
\name{result.DEG}
\alias{result.DEG}
\title{Get a summary result from multiple DEGList}
\description{
result.DEG help get a summary result from multiple DEGList
}
\usage{
rd <- result.DEG(DEGList,
                 logFC.list,
                 select,
                 filter=T)
}
\arguments{
\item{DEGList}{a list of result from DESeq1 or edgeR1}
\item{logFC.list}{a list of logFC colnames}
\item{select}{selected markers or genes}
\item{filter}{whether to do a same trend filter.If you want all results,set FALSE to filter}
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
## This is a simulative process and NOT RUN

# parameters
DEGList <- list(DESeq2 = train.DESeqList,edgeR = train.edgeRList)
logFC.list <- list(DESeq2 = "log2FoldChange",edgeR = "logFC")
select <-  intersect(train.DESeqList$siggenes,train.edgeRList$siggenes)

# result
resultDEG <- result.DEG(DEGList,
                        logFC.list,
                        select,
                        filter = T)


}
\keyword{result.DEG}




