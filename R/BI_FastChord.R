

FastChord <- function(data,
                      source=NULL,
                      target=NULL,
                      value=NULL,
                      grid.col = NULL,
                      transparency = 0.5,
                      col = NULL,
                      save.path = "Chord",
                      names="love"){

  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)

  ## 区分邻接矩阵(adjacency matrix)和邻接列表(adjacency list)
  if(is.null(source)|is.null(target)|is.null(value)){
    LuckyVerbose("data is an adjacency matrix...")
    data_1 <- data
  } else {
    LuckyVerbose("data is an adjacency list...")
    data_1 <- data.frame(
      from = as.character(data[,source]),
      to = as.character(data[,target]),
      value = as.numeric(as.character(data[,value]))
    )
  }

  ## draw chord diagram
  #win.graph(12,12)
  chordDiagram(x = data_1 ,
               col = col,
               grid.col = grid.col,
               transparency = transparency)

}















