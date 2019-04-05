
\name{FastMINE}
\alias{FastMINE}
\title{a fast way to do MINE correlation analysis.}
\description{
 FastMINE通过MINE策略挖掘变量间的相关性。MINE不需要考虑变量的分布，在global health/gene expression/major-league baseball/the human gut microbiota等大数据挖掘面前有一定的优势。不过当背景噪音很强时，MINE的表现出相对较弱的相关性挖掘能力。但总体上，在相关性挖掘方面，MINE算法要优于很多传统的算法(比如Pearson法)。FastMINE()的运行依赖Java。
}
\usage{
  FastMINE(data,
         transposition = F,#是否转置矩阵
         control.markers,
         target.markers=NULL,
         method=c("one.pair","all.pairs")[2])
}
\arguments{

}
\details{
某些可以考虑MINE策略的特殊情况：
noncoexistence relationships(非共存关系):一者存在，另一者几乎不存在。


}
\value{
\item{X var/Y var}{作比较的markers pair.}
  \item{MIC (strength)}{the maximal information coefficient.MIC is general and roughly equal to R2 on functional relationships.相关性强度。MIC越高，相关性越强。与Pearson correlation是类似的。}
  \item{MIC-p^2 (nonlinearity)}{MIC-Pearson^2，非线性程度。y=x的MIC-p^2为0。MIC-p^2值越高，相关性的非线性程度就越强。}
  \item{MAS (non-monotonicity)}{非单调性。y=x的MAS为0。MAS越高，则非单调性越强,越难用单调函数进行相关性描述。}
  \item{MEV (functionality)}{函数性。y=x的MEV为1,MEV越高代表函数性越强，越容易用函数进行描述。}
  \item{MCN (complexity)}{复杂度:y=x的值为2，越高表明相关性越复杂，越难用尽量少的变量去描述它们之间的相关性。}
  \item{Pearson(p)}{Pearson correlations}
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
library(lucky)
data=state.x77;colnames(data)
result1 <- FastMINE(data,
                    transposition = F,
                    method = "all.pairs",
                    control.markers="Income",
                    target.markers=NULL)
View(result1[["MIC.matirx"]])
result2 <- FastMINE(data,
                    transposition = F,
                    method = "one.pair",
                    control.markers="Income",
                    target.markers=NULL)
colnames(result2);View(result2)
}
\keyword{FastMINE}




