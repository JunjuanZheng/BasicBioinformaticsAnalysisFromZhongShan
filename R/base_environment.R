

####========================lucky全局变量=========================####

###======自动加载需要的包
# package.need <- c("stringr","ggplot2")
# for(i in package.need){require(i,character.only = T)}

###======mymusic

#' @title play music via tuneR
#' @description play music via tuneR
#' @keywords mymusic
#' @param n the order of available music paths
#' @param path the path of music files ducument
#' @return a string of music file path
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' mymusic(1)
#' @export
mymusic <- function(n=1,path=NULL){
  if(is.null(path)){
    path <- "E:/RCloud/music"
  }
  music.list <- list.files(path,pattern = "Rend_",full.names = T)
  if(length(music.list)==0){
    # 在给定的地址中未找到音频文件。
    library(tuneR)
    music <- tuneR::bind(sine(440), sine(220))
    tuneR::play(music)
    #print("未找到音频，自动生成一个简单音频。")
  } else {
    # 可以找到文件
    #print(paste0("一共有",length(music.list),"个音频文件。"))
    music <- music.list[n];
    tuneR::play(music)
    #print(paste0("选择了第",select,"个音频文件。"))
  }
}






