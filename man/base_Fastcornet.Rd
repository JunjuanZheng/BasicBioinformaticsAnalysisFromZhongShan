
\name{Fastcornet}
\alias{Fastcornet}
\title{fast way to create a network plot and edge/node data based on a cor matrix.}
\description{
 fast way to create a network plot and edge/node data based on a cor matrix.
}
\usage{
Fastcornet(cormatrix,
           color.group,
           ratio,
           legend=T,
           size = 15,
           control.marker=NULL,
           control.linkColour="red",
           other.linkColour=NULL,
           linkColour.type = list(
             name = "YlGn",
             select.lower = 1,
             select.upper = 9
           ),
           colourScale.type = 1,
           colourScale.index = 1,
           linkDistance = 300,
           linkWidth=6,
           arrows=F,
           zoom = T)
}
\arguments{
  \item{cormatrix}{a matrix of correlation index.}
  \item{color.group}{a named list of markers categery}
  \item{control.group}{the names of control group}
  \item{ratio}{select the most significant portion of correlations.it range 0 to 1.}
  \item{legend}{whether show legend in the netplot output}
  \item{size}{the index to enhance the degree difference among nodes.If you consider the size of nodes too small,try make size larger.}
  \item{linkColour}{link colour}
  \item{linkColour.type}{a list of parameters to RColorBrewer::brewer.pal}
  \item{control.linkColour}{the color of links to control group nodes.If NULL,use default colours,which is recommanded}
  \item{colourScale.type}{d3-scale-chromatic style object.1="d3.schemeCategory20",2="d3.schemeCategory20b",3="d3.schemeCategory20c",4="d3.schemeCategory10".Default is 1}
  \item{nodeColour}{the named list of group color.It is unavailable when colourScale.type exists}
  \item{linkDistance}{link distance}
  \item{linkWidth}{the width of link.If it is not a numeric,the default "JS(function(d) { return Math.sqrt(d.value); })" would be used.}
  \item{opacity}{numeric value of the proportion opaque you would like the graph elements to be}
  \item{arrows}{F.whether show direction between two nodes}
  \item{zoom}{logical value to enable (TRUE) or disable (FALSE) zooming}

}
\details{
  linkColour.type make full use of RColorBrewer and would be a good color strategy provider.The parameters of linkColour.type is including:
  1.names: One of "Blues","BuGn","BuPu","GnBu","Greens","Greys","Oranges","OrRd","PuBu","PuBuGn","PuRd","Purples","RdPu","Reds" ,"YlGn","YlGnBu","YlOrBr" and "YlOrRd";
  2.select.lower: the lower color order;
  3.select.upper: the upper color order.
}
\value{
 a list contain a network plot information and edge/node data.
}
\references{
 networkD3 package
}
\author{
  Weibin Huang<654751191@qq.com>
}
\note{

}
\seealso{
 networkD3::forceNetwork
 lucky::D3.colourScale
}
\examples{
## example
library(lucky)
data=state.x77;colnames(data)
## get MIC cormatrix
result1 <- FastMINE(data,
                    transposition = F,
                    method = "all.pairs",
                    control.markers="Income",
                    target.markers=NULL)
cormatrix <- result1[["MIC.matirx"]]

## parameters
# group
color.group = list(
  love = c("Murder","Frost","Area","HS Grad"),
  hate = c("Life Exp","Illiteracy"),
  control = "Population")
# group color
nodeColour=list(
  love = mycolor[1],
  hate = mycolor[5],
  control = mycolor[4])
# link color
linkColour.type = list(
  name = "YlGn",
  select.lower = 1,
  select.upper = 9
)

# control.group = "control";control.linkColour = "red"
control.group = NULL;control.linkColour =NULL
colourScale.type = NULL
# colourScale.type = 1

## Quick Start
nw1 <- Fastcornet(cormatrix,
                  color.group,
                  control.group,
                  ratio=0.8,
                  legend=T,
                  size = 15,
                  linkColour=NULL,
                  linkColour.type,
                  control.linkColour,
                  colourScale.type,
                  nodeColour);nw1$netplot

}
\keyword{Fastcornet}

