
\name{result.lassoModel}
\alias{result.lassoModel}
\title{Get results from lasso model function}
\description{
 Get results from lasso model function
}
\usage{
result.lassoModel(lassoModel,
                  position = c(2,1,1),
                  names = c("TTR","DFS","OS"),
                  dig = 5)
getModelMatrix(expr.matrix,
               design,
               model,
               cut.off = NULL,
               position = 1)

}
\arguments{

1.result.lassoModel:
  \item{lassoModel}{the result from cox.optimized2 and bootstrap.lasso series}
  \item{position}{the position of one model you want}
  \item{names}{the name of one model like DFS}
  \item{dig}{digits of correlation of model genes}

2.getModelMatrix:
  \item{expr.matrix}{a log scale expression matrix}
  \item{design}{a design object}
  \item{model}{result from result.lassoModel}
  \item{cut.off}{Default is NULL.If you set a numeric,it would produce a score status column}
  \item{position}{the position of one model you want}
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

model <- result.lassoModel(l1,
                           position = 1,
                           names = "ROC",
                           dig = 5)
View(model[[1]])

## get model matrix
model.mt <- getModelMatrix(expr.matrix = data.fpkm,
                           design = design.train,
                           model = model)
View(model.mt$metadata)
}
\keyword{result.lassoModel}




