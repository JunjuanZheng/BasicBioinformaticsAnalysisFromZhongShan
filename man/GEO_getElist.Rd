
\name{getElist}
\alias{getElist}
\title{getElist}
\description{
 getElist()通过给定GSE号，自动下载标准格式的soft文件，提取其中的矩阵、注释和target信息，最后形成eset类列表。此函数依赖于GEOquery，可以获得详尽的表型、注释、处理流程内容，十分强大。不过，如果是外地下载soft文件，应该改为GSEXXXXX.soft.gz的名称才能正确识别。如果想保存eset文件为rda，可以选择savefile = T。如果有个性化GSE.soft文件保存空间，可以自定义gse.space属性。
}
\usage{
  getElist(gse.names,
         gse.space=NULL,
         savefile = F)
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

}
\keyword{compareROC2}

