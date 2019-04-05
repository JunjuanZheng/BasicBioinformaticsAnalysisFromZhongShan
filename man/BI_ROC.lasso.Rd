
\name{ROC.lasso}
\alias{ROC.lasso}
\title{Bulid a binomial lasso model via exprs and design objects}
\description{
ROC.lasso help bulid a binomial lasso model via exprs and design objects.
}
\usage{
ROC.lasso(expr.matrix,
          design,
          select,
          contrast = "N.status",
          contrast.control = "N0",
          k=10,
          R=100,
          seed = 2018,seed.range = 2500:5000,
          optimize.method = "min",
          show.music = F)
}
\arguments{
  \item{expr.matrix}{a log scale expression matrix}
  \item{design}{a design object}
  \item{select}{select genes names}
  \item{contrast}{the colnames of the contrast}
  \item{contrast.control}{the control of the contrast}
  \item{k}{nfolds}\\item{R}{the number of bootstrap}
  \item{seed}{random seed}
  \item{seed.range}{the range of random seeds}
  \item{optimize.method}{one of "min" and "1se".Default is "min"}
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

library(lucky)
load1(c("data.fpkm","design.train"))

## main parameters
expr.matrix = data.fpkm
design = design.train
select = rownames(data.fpkm)[1:50]

## others parameters
contrast = "N.status"
contrast.control = "N0"
k=10
R=100
seed = 2018
seed.range = 2500:5000
each.size = 15
optimize.method = "min"
show.music = T

## ROC.lasso
l1 <- ROC.lasso(expr.matrix,
                design,
                select,
                contrast = "N.status",
                contrast.control = "N0",
                k=10,
                R=100,
                seed = 2018,seed.range = 2500:5000,
                optimize.method = "min",
                show.music = F)

## look at the models
View(uniqueModel(l1$modeldata))

## select a model in position 1
model <- result.lassoModel(l1,
                           position = 1,
                           names = "ROC",
                           dig = 5)
View(model[[1]])

## get model matrix for ROC building
x <- getModelMatrix(expr.matrix = data.fpkm,
                    design = design.train,
                    model = model)
View(x$metadata)

## ROC building
roc1 <- compareROC2(data = x$metadata,
                    markers = c(x$modelgenes,"score"),
                    status = "N.status",
                    merge.markers=list(merge = x$modelgenes),
                    roc.type = "ggplot",
                    title="ROC Curves",
                    color=NULL,
                    half.border = F,reference = T,
                    legend.position=c(0.8,0.2),
                    show.auc = T,auc.digits = 2,
                    width = "love")
}
\keyword{ROC.lasso}




