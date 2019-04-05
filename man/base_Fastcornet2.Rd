
\name{Fastcornet2}
\alias{Fastcornet2}
\title{fast way to create a network plot and edge/node data based on a long data frame from a cormatrix.}
\description{
 fast way to create a network plot and edge/node data based on a long data frame from a cormatrix.
}
\usage{
Fastcornet2(data,
            source.col=NULL,target.col=NULL,value=NULL,
            color.group = list(),
            control.group = "control",
            ratio=0.8,
            legend=T,
            size = 15,
            linkColour=NULL,
            linkColour.type = list(
              name = "YlGn",
              select.lower = 1,
              select.upper = 9
            ),
            control.linkColour = "red",
            colourScale.type = NULL,
            nodeColour=list(),
            linkDistance = 300,
            linkWidth=6,
            opacity=1,
            arrows=F,
            zoom = T)
}
\arguments{
  \item{data}{a long data frame from a cormatrix}
  \item{source.col,target.col,value}{parameters as Make.cormatrix() in lucky package.}
  \item{color.group}{a named list of markers categery}
  \item{control.group}{the names of control group}
  \item{ratio}{select the most significant portion of correlations.it range 0 to 1.}
  \item{legend}{whether show legend in the netplot output}
  \item{size}{the index to enhance the degree difference among nodes.If you consider the size of nodes too small,try make size larger.}
  \item{linkColour}{link colour}
  \item{linkColour.type}{a list of parameters to RColorBrewer::brewer.pal}
  \item{control.linkColour}{the color of links to control group nodes.If NULL,use default colours,which is recommanded}
  \item{colourScale.type}{d3-scale-chromatic style object.1="d3.schemeCategory20",2="d3.schemeCategory20b",3="d3.schemeCategory20c",4="d3.schemeCategory10".Default is 1}
  \item{nodeColour}{}
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
 networkD3::forceNetwork()
 lucky::Fastcornet
}
\examples{
library(lucky)
data=state.x77;colnames(data)
## get MIC cormatrix
result1 <- FastMINE(data,
                    transposition = F,
                    method = "all.pairs",
                    control.markers="Income",
                    target.markers=NULL)

## long data frame network
data <- result1[["MINE.result"]][,1:3]
colnames(data) #[1] "X var"  "Y var" "MIC (strength)"
source.col="X var";target.col="Y var";value="MIC (strength)"

## parameters of Fastcornet2
color.group = list(
  love = c("Murder","Frost","Area","HS Grad"),
  hate = c("Life Exp","Illiteracy"),
  control = c("Population","Income"))
control.group = "control"
linkColour.type = list(
  name =  "Blues",
  select.lower = 1,
  select.upper = 3
)

# control.linkColour = "red"
control.linkColour = NULL
colourScale.type = NULL
nodeColour=list(
  love = mycolor[6],
  hate = mycolor[8],
  control = mycolor[35])

## Quick Start
nw2 <- Fastcornet2(data,
                   source.col,
                   target.col,
                   value,
                   color.group,
                   control.group,
                   ratio=0.8,
                   legend=T,
                   size = 15,
                   linkColour=NULL,
                   linkColour.type,
                   control.linkColour,
                   colourScale.type = NULL,
                   nodeColour);nw2$netplot
}
\keyword{Fastcornet2}

