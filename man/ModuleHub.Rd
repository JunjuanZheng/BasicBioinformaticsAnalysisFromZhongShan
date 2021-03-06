% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BI_ModuleHub.R
\name{ModuleHub}
\alias{ModuleHub}
\title{Enhanced hub genes exploration after FastWGCNA pipeline}
\usage{
ModuleHub(object, design, variable, corType = "pearson",
  cutoff.pval = 1, hub_cutoffSigGM = 0.2, hub_MM = 0.8,
  hub_WeightedQ = 1, save.path = "WGCNA", names = "love")
}
\arguments{
\item{object}{the result of  \code{\link{FastWGCNA}}.}

\item{design}{a trait-design object}

\item{variable}{the variables you want to show in Module-Trait relationships plot.}

\item{corType}{one of "pearson" and "bicor".Default is pearson}

\item{cutoff.pval}{cut-off of the p value in significant module-phenotype filter}

\item{hub_cutoffSigGM}{the cut-off of significant genes-Modules in hub genes exploration.Default is 0.2}

\item{hub_MM}{the cut-off of Module Memberships in hub genes exploration.Default is 0.8}

\item{hub_WeightedQ}{the cut-off of weighted q value in hub genes exploration.Default is 0.01}

\item{save.path}{the space of the save file.Default is "WGCNA"}

\item{names}{part of saved files name}
}
\value{
LuckWGCNA object
}
\description{
The speed of \code{ModuleHub} is very fast.If you want to set different parameters(always \code{hub_WeightedQ} and \code{cutoff.pval}) for hub genes exploration, \code{ModuleHub} is a nice function to do it.
}
\details{
1.Only work on wgcna result with ONE Block.   2.\code{cutoff.pval}
  is a useful parameter.If you want significant module,you can set
  \code{cutoff.pval = 0.05}(or 0.01,It depends on your custom);If you want to
  see all the modules regardless of significance,just set \code{cutoff.pval =
  1}.  3.\code{hub_WeightedQ} is a stricter filter for hub genes.If your hub
  genes is too much,you can set \code{hub_WeightedQ = 0.05}(or 0.01,It
  depends on your custom).However,most researchers do not use
  \code{hub_WeightedQ} to filter their hub genes and often use only
  \code{hub_MM=0.8} and \code{hub_cutoffSigGM=0.2}.
}
\examples{
## This is a simulative process and available only with CORRECT VARIABLES
library(lucky)
load("E:/RCloud/RFactory/lucky/love/WGCNA-test/love_wgcna.rda")
object = wgcna;rm(wgcna);gc()
design = rna.design.tumor
variable = c("age","his1","gender","N.status","T.status")
result_MH <- ModuleHub(object,
                       design,
                       variable,
                       save.path = "WGCNA",
                       names = "love")
}
\seealso{
\code{\link{FastWGCNA}}.
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
\keyword{ModuleHub}
