
\name{bootstrap.lasso}
\alias{bootstrap.lasso}
\title{bulid a bootstrap lasso model via exprs and design objects}
\description{
bulid a bootstrap lasso model for binary induced variable via exprs and design objects and output a ROCList object
}
\usage{
bootstrap.lasso(dds,
                transformation = "normTransform",
                expr.matrix=NULL,
                design,
                select,
                contrast = "N.status",
                contrast.control = "N0",
                k=10,
                R=2000,
                seed = 2018,seed.range = 1:15000,
                each.size = 15,
                optimize.method = "1se",
                show.music = T)

View(uniqueModel(lassoModel$modeldata))

result.lassoModel(lassoModel,
                  position = 1,
                  names = "test",
                  dig = 5)
}
\arguments{
  \item{dds}{DESeq2 dds object.If it was NULL,the expr.matrix should not be null.}
  \item{transformation}{one of "vst" and "normTransform".}
  \item{expr.matrix}{a log scale expression matrix.}
  \item{design}{a design object.}
  \item{select}{select genes names.}
  \item{contrast}{the colnames of the contrast}
  \item{contrast.control}{the control of the contrast}
  \item{k}{nfolds}
  \item{R}{2000.bootstrap}
  \item{seed}{random seed}
  \item{seed.range}{the range of random seeds}
  \item{each.size}{the size of each group in one bootstrap}
  \item{optimize.method}{one of "min" and "1se".Default is "1se".}
  \item{show.music}{whether show music at the end of the function}
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
result.lassoModel() function also support the result from bootstrap.lasso
}
\seealso{

}
\examples{
## a virtual not-running programes
ROCList <- bootstrap.lasso(expr.matrix=model1.fpkm,
                           design = design.model1,
                           select = LM.genes,
                           contrast = "N.status",
                           contrast.control = "N0",
                           k=10,
                           R=2000,
                           seed = 2018,seed.range = 1:15000,
                           each.size = 15,
                           optimize.method = "1se")
}
\keyword{bootstrap.lasso}




