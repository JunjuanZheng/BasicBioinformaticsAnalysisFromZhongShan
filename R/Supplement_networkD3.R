



# id # type is one of Categorical,Diverging,Sequential.mh and Sequential.sh and the order is an integer.like Categorical-1
# report # whether to show all available colour scales
#' @export
D3.colourScale <- function(id = "Categorical-1",
                           report=F){
  ### All color
  {
  ## Categorical
  Categorical <- c("d3.schemeCategory20",
                   "d3.schemeCategory10",
                   "d3.schemeAccent",
                   "d3.schemeDark2",
                   "d3.schemePaired",
                   "d3.schemePastel1",
                   "d3.schemePastel2",
                   "d3.schemeSet1",
                   "d3.schemeSet2",
                   "d3.schemeSet3")
  names(Categorical) <- paste("Categorical",1:length(Categorical),sep = "-")

  ## Diverging
  Diverging <- c("d3.interpolateBrBG",
                 "d3.interpolatePiYG",
                 "d3.interpolatePRGn",
                 "d3.interpolatePuOr",
                 "d3.interpolateRdBu",
                 "d3.interpolateRdGy",
                 "d3.interpolateRdYlBu",
                 "d3.interpolateRdYlGn",
                 "d3.interpolateSpectral",
                 "d3.schemeBrBG",
                 "d3.schemePiYG",
                 "d3.schemePRGn",
                 "d3.schemePuOr",
                 "d3.schemeRdBu",
                 "d3.schemeRdGy",
                 "d3.schemeRdYlBu",
                 "d3.schemeRdYlGn",
                 "d3.schemeSpectral")
  names(Diverging) <- paste("Diverging",1:length(Diverging),sep = "-")

  ## Sequential (Single Hue)
  Sequential.sh <- c("d3.interpolateBlues",
                  "d3.interpolateGreens",
                  "d3.interpolateGreys",
                  "d3.interpolateOranges",
                  "d3.interpolatePurples",
                  "d3.interpolateReds",
                  "d3.schemeBlues",
                  "d3.schemeGreens",
                  "d3.schemeGreys",
                  "d3.schemeOranges",
                  "d3.schemePurples",
                  "d3.schemeReds")
  names(Sequential.sh) <- paste("Sequential.sh",1:length(Sequential.sh),sep = "-")

  ## Sequential (Multi-Hue)
  Sequential.mh <- c("d3.interpolateBuGn",
                     "d3.interpolateBuPu",
                     "d3.interpolateCool",
                     "d3.interpolateCubehelixDefault",
                     "d3.interpolateGnBu",
                     "d3.interpolateInferno",
                     "d3.interpolateMagma",
                     "d3.interpolateOrRd",
                     "d3.interpolatePlasma",
                     "d3.interpolatePuBu",
                     "d3.interpolatePuBuGn",
                     "d3.interpolatePuRd",
                     "d3.interpolateRdPu",
                     "d3.interpolateViridis",
                     "d3.interpolateWarm",
                     "d3.interpolateYlGn",
                     "d3.interpolateYlGnBu",
                     "d3.interpolateYlOrBr",
                     "d3.interpolateYlOrRd",
                     "d3.schemeBuGn",
                     "d3.schemeBuPu",
                     "d3.schemeGnBu",
                     "d3.schemeOrRd",
                     "d3.schemePuBu",
                     "d3.schemePuBuGn",
                     "d3.schemePuRd",
                     "d3.schemeRdPu",
                     "d3.schemeYlGn",
                     "d3.schemeYlGnBu",
                     "d3.schemeYlOrBr",
                     "d3.schemeYlOrRd")
  names(Sequential.mh) <- paste("Sequential.mh",1:length(Sequential.mh),sep = "-")
}
  all <- c(Categorical,Diverging,Sequential.mh,Sequential.sh)

  ### Report
  if(report == T){
    print(paste0("&====================1.Categorical: ",length(Categorical)," type===============&"))
    print(Categorical)
    print(paste0("&====================2.Diverging: ",length(Diverging)," type=================&"))
    print(Diverging)
    print(paste0("&==================3.Sequential.sh: ",length(Sequential.sh)," type===============&"))
    print(Sequential.sh)
    print(paste0("&==================4.Sequential.mh: ",length(Sequential.mh)," type===============&"))
    print(Sequential.mh)
  }

  ### Select an object
  lg1 <- names(all) %in% id
  if(T %in% lg1){
    ## 说明存在相应编号
    one <- all[lg1]
  } else {
    one <- "Error:Not such D3.colourScale ID!"
  }

  return(as.character(one))

}






