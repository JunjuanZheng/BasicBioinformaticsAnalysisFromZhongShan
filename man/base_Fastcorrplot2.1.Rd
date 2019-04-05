
\name{Fastcorrplot2.1}
\alias{Fastcorrplot2.1}
\title{correlation analysis limited one control marker}
\description{
 Give one control markers then test the correlation between every control marker and target markers.Output a corrplot based on lucky::Fastcorrplot() function and a data frame containing compaired pairs and P values.
}
\usage{
  Fastcorrplot2.1(data,
                transposition = T,#是否转置矩阵
                control.markers,
                target.markers=NULL,
                method="pearson",p.cut.off=0.05,
                savefile=T,#corrplot()相关参数
                names="test1",
                lower.col = NULL,#corrplot()相关参数
                upper.col =NULL,#corrplot()相关参数
                upper = NULL,#corrplot()相关参数
                tl.pos = NULL,#corrplot()相关参数
                tl.col=NULL,#corrplot()相关参数
                tl.srt=NULL,#corrplot()相关参数
                diag = NULL)
}
\arguments{
  \item{data}{a gene expression matrix or a data frame with patient cols and gene rows.Or a matrix and a dataframe with similar construction.}
  \item{transposition }{ T.whethe make data transposition.}
  \item{control.markers}{a gene used as internal reference to explore correlations.}
  \item{target.markers}{NULL.other genes differing from control markers.It must be part of data rownames(transposition=T) or colnames.If NULL,it is the other set beyond control markers.}
  \item{method}{please see psych::corr.test}
  \item{p.cut.off}{a cut.off of significance.Default=0.05.}
  \item{savefile}{T.Whether to save a PDF plot.}
  \item{names}{part of PDF file name.}
  \item{other parameters}{please see corrplot::corrplot().}
}
\details{
    method参数来自psych::corr.test函数中的method参数。选择条件如下：
    1.两个连续变量间呈线性相关时，使用Pearson积差相关系数，不满足积差相关分析的适用条件时，使用Spearman秩相关系数来描述.
    2.Spearman相关系数又称秩相关系数，是利用两变量的秩次大小作线性相关分析，对原始变量的分布不作要求，属于非参数统计方法，适用范围要广些。对于服从Pearson相关系数的数据亦可计算Spearman相关系数，但统计效能要低一些。Pearson相关系数的计算公式可以完全套用Spearman相关系数计算公式，但公式中的x和y用相应的秩次代替即可。
    3.Kendall's tau-b等级相关系数：用于反映分类变量相关性的指标，适用于两个分类变量均为有序分类的情况。对相关的有序变量进行非参数相关检验；取值范围在-1-1之间，此检验适合于正方形表格
    如果不知道选哪个，就选择默认的Pearson法即可。
}
\value{
 a corrplot + a data frame
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
## example
data=state.x77
control.markers="Income";
target.markers=NULL
names = "state.x77"
cor.df <- Fastcorrplot2.1(data,
                          transposition=F,#是否转置矩阵
                          control.markers,
                          target.markers=NULL,
                          p.cut.off=0.05,
                          savefile=T,#corrplot()相关参数
                          names=names)

## change cor calculated method
cor.df <- Fastcorrplot2.1(data,
                          transposition=F,#是否转置矩阵
                          method = "spearman",
                          control.markers,
                          target.markers=NULL,
                          p.cut.off=0.05,
                          savefile=T,#corrplot()相关参数
                          names=names)
}
\keyword{Fastcorrplot2.1}




