




####====Development tools for Weibin Huang===####
####=======创建依赖包的关系
#' @export
make.dependent.packages <- function(packages,type=1){
  nd <- c("devtools","stringr");Plus.library(nd)
  packages <- sort(packages)
  types=c("Suggests","Imports")
  for(i in packages){
    use_package(i,types[type])
  }
  # End
}
# sug <- c("rJava","xlsx","dplyr") type=1即Suggests
# imp <- c("stringr","plyr")# type=2即"Imports"
# make.dependent.packages(imp,1);make.dependent.packages(sug,2)


####=======创建Rb文件中的arguement
## 作用：将source文件中的备注转化为Rd文件中的arguments.
# select.genes=a #待处理的基因。
# design=b #含有分组信息的数据框
#' @export
make.arguments <- function(remark){
  nd <- c("devtools","stringr");Plus.library(nd)
  r1 <- unlist(strsplit(remark,"#"))
  ## 将换行符去除。
  enter.symbol.c <- grep("\n",r1)
  r2 <- r1[enter.symbol.c]
  r2.1 <- NULL
  for(i in r2){
    r2.1i <- unlist(strsplit(i,"\n"))
    r2.1 <- c(r2.1,r2.1i)
  }
  r1[enter.symbol.c] <- r2.1

  ## 将空字符去除
  r1.nchar <- nchar(r1)
  r3 <- r1[r1.nchar != 0]

  ## 将字符首尾的空格键删除
  #vt=r3
  #vt=rep("merge.i.values=list(merge1 = i.values) ",3)
  #end.start.space.delet(vt)
  end.start.space.delet <- function(vt){
    if(length(vt)==1){
      #仅一个元素
      # unlist(strsplit(" vt "," "))
      vt1 <- unlist(strsplit(vt," "))
      vt1 <- vt1[nchar(vt1) !=0]
    } else {
      vt1 <- NULL
      #i = " color = brewer.pal(12, \"Set3\")[1:5] "
      #i=" merge.markers = list(c(\"BMI\",\"Waist\"))"
      #i= "merge.i.values=list(merge1 = i.values) "
      for(i in vt){
        vt.i <- unlist(strsplit(i," "))
        vt.i <- vt.i[nchar(vt.i) !=0]
        vt.i <- paste0(vt.i,collapse = " ")
        vt1 <- c(vt1,vt.i)
      }
    }
    return(vt1)
  }#去除字符首尾的空字符
  r4.1 <- end.start.space.delet(r3)

  ## 构建arguement类字符串
  #str.i = c("abc","bcd")
  #str.i = c("a","b")
  #str.i = c("a=list(a=b)","a list")
  argue.i <- function(str.i){
    l1 <- grep("=",str.i);l1 <- length(l1)
    if(l1 == 0){
      ## 说明没有等号。可以直接构建argument.
      str.ix1 <- str_c("\\","item","{",str.i[1],"}","{",str.i[2],"}")
    } else {
      ## 说明有等号。判断某个字符里有几个=。只删除第一个=
      #str="a=list(a=b)"
      #str="a=list(ab)"
      #str="a"
      #str="a list"
      del.1 <- function(str){
        str.m <- unlist(strsplit(str,"="))
        l1 <- length(str.m)-1
        if(l1==0){
          #没有等号
          str1 <- str
        } else {
          if(l1==1){
            #1个等号
            str1 <- str.m
          } else {
            #2个及以上的等号。只去除第一个等号
            str1 <- c(str.m[1],paste0(str.m[2:length(str.m)],collapse = "="))
          }
        }
        return(str1)
      }#只删除第一个=,生成2个字符del.1(str)
      str.ix <- NULL
      for(i in str.i){
        #i = "a list"
        str.ix1 <- del.1(i)
        str.ix <- c(str.ix,str.ix1)
      }
      str.ix <- c(str.ix[1],paste0(str.ix[2:length(str.ix)],collapse = "."))
      str.ix1 <- str_c("\\","item","{",str.ix[1],"}","{",str.ix[2],"}")
    }
    return(str.ix1)
    # median end.
  } #argue.i(str.i)
  #i=1
  argue <- NULL
  for(i in 1:(length(r4.1)/2)){
    str.i <- c(r4.1[2*i-1],r4.1[2*i])
    a.i <- argue.i(str.i)
    argue <- str_c(argue,a.i)
  }
  argue

  ##输出结果
  return(argue)
}
## 例：备注规范：每行只有两个"#"字符。=号两边没空格。
# remark <- "# select.genes=a #待处理的基因。
# design=b #含有分组信息的数据框"
# make.arguments(remark)
# 输出："\\item{select.genes}{a.待处理的基因。}\\item{design}{b.含有分组信息的数据框}"




