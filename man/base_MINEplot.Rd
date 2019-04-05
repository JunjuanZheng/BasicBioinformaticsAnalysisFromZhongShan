
\name{MINE.corrplot}
\alias{MINE.corrplot}
\title{a corrplot based on MINE cor matrix}
\description{
 a corrplot based on MINE cor matrix
}
\usage{
  MINE.corrplot(data,
              transposition = F,
              select,
              order = T,
              lower.col = NULL,
              upper.col =NULL,
              upper = NULL,
              tl.pos = NULL,
              tl.col=NULL,
              tl.srt=NULL,
              diag = NULL,
              savefile=T,
              names="love")
}
\arguments{
  \item{data}{gene expression matrix or data frame with patient cols and gene rows.}
  \item{transposition}{Whether use t(data) before corrplot making.If your matrix is one with gene rows,you should let transposition=T.}
   \item{order}{Whether use AOE order strategy.Default is T.}
  \item{select}{genes your want to explore the correlations.Note that the select genes must be part of rownames of expr.}
  \item{...}{parameters from corrplot.mixed.Default = NULL.Also you can defined it alternatively.}
  \item{savefile}{whether saving plot as PDF file.}
  \item{names}{part of PDF file names.}
}
\details{

}
\value{
 corrplot
 a list containing cor matrix of MIC values.
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
## example
data = state.x77
result <- MINEplot(data,
                   transposition = F,
                   select=colnames(data),
                   QC=T,order = T,outer.ratio = 0.5,
                   savefile=F,
                   names="love")
View(result[["MINE.result"]])
View(result[["MIC.matirx"]])

}
\keyword{MINE.corrplot}




