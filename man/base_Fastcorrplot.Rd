
\name{Fastcorrplot}
\alias{Fastcorrplot}
\title{help fast to draw a corrplot with multiple styles available to select}
\description{
 help fast to draw a corrplot with multiple styles available to select
}
\usage{
 Fastcorrplot(data,
             transposition = T,#是否转置矩阵
             order = T,
             select,
             lower.col = NULL,
             upper.col =NULL,
             upper = NULL,
             tl.pos = NULL,
             tl.col=NULL,
             tl.srt=NULL,
             diag = NULL,
             width=10,height=10,
             savefile=T,
             names="test1")
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
Fastcorrplot() help fast to draw a corrplot with multiple styles available to select.Fastcorrplot() is based on corrplot::corrplot.mixed().
}
\value{
 a cor matrix and a corrplot
}
\references{
corrplot.mixed()
}
\author{
  Weibin Huang
}
\note{

}
\seealso{
corrplot.mixed()
}
\examples{
data(expr)
set.seed(2018);select=rownames(expr)[sample(1:nrow(expr),10)]
data <- expr[,1:(ncol(expr)-4)]
cor.matrix <- Fastcorrplot(data,
                           transposition = T,
                           select,
                           savefile=F,
                           names="test1")
# self-defined upper parameter
upper <- c("ellipse","circle","square","shade","pie")
cor.matrix <- Fastcorrplot(data,
                           transposition = T,
                           upper = upper[1],
                           select,
                           savefile=F,
                           names="test1")

# self-defined transposition parameter
upper <- c("ellipse","circle","square","shade","pie")
cor.matrix <- Fastcorrplot(data=state.x77,
                           transposition = F,
                           upper = upper[5],
                           select = colnames(state.x77),
                           savefile=F,
                           names="test1")
}
\keyword{Fastcorrplot}




