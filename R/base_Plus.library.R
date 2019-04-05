

####======================Plus.library()=======================####
#自动下载、安装、加载包。
# packages = c("ggplot2","MASS","gplots")

#'Automatically install and library package
#'@description A fast way for Chinese people to automatically install and library package.
#'@param packages a string vector of package names
#'
#'@author Weibin Huang\email{654751191@@qq.com}
#'
#'@examples
#'Plus.library(c("DT","RColorBrewer","stringr"))
#'
#'@export
Plus.library <- function(packages){
  ## 系统已经安装的包
  x <- installed.packages()
  x <- as.data.frame(x)
  installed <- as.character(x$Package)

  ## 检查系统中有没有安装pacman包
  is <- length(grep("pacman",installed))
  if(is ==0){
    #pacman包未安装
    install.packages("pacman")
    library(pacman)
  } else {
    library(pacman)
  }

  ## 加载包
  for(i in packages){p_load(char=i)}
}


