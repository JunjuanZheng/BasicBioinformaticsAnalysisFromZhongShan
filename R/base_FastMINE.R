


###=====FastMINE：a fast way to do MINE correlation analysis.
## Usage FastMINE通过MINE策略挖掘变量间的相关性。MINE不需要考虑变量的分布，在global health/gene expression/major-league baseball/the human gut microbiota等大数据挖掘面前有一定的优势。不过当背景噪音很强时，MINE的表现出相对较弱的相关性挖掘能力。但总体上，在相关性挖掘方面，MINE算法要优于很多传统的算法(比如Pearson法)。FastMINE()的运行依赖Java。
#' @export
FastMINE <- function(data,
                     transposition = F,#是否转置矩阵
                     control.markers,
                     target.markers=NULL,
                     method=c("one.pair","all.pairs")[1]){
  ## 加载包
  # print("Fastcorrplot.MINE depends on the correct installation of JAVA in your computer.If Fastcorrplot.MINE not work,please go online for JAVA installation in Windows platform.")
  library(rJava);

  ### data transposition
  expr1 <- as.matrix(data)
  if(transposition==T){
    # 矩阵要转置
    expr2 <- t(expr1)
  } else {
    expr2 <- expr1
  }

  ## 加载默认设置
  l1 <- colnames(expr2) %in% control.markers
  default <- c(rep(list(NULL),1)) #
  input <- list(target.markers=target.markers)
  do <- list(target.markers=colnames(expr2)[!l1])
  output <- set.default(input,default,do);
  expr2 <- expr2[,c(control.markers,output[[1]])]
  #names(output)
  #[1] "target.markers" "lower.col"      "upper.col"
  #[4] "upper"          "tl.pos"         "tl.col"
  #[7] "tl.srt"         "diag"

  ### 将目标数据保存到一个现实csv文件
  csv.path <- "./MINE.test.csv"
  print("data.table:正在本地写入csv文件...")
  library(data.table)
  dt.expr <- as.data.frame(expr2)
  fwrite(dt.expr,file = csv.path,row.names = F)
  print("data.table:完成本地写入csv文件!")

  ### 加载java VM的位置.依赖于Supplement_MINE.R
  java.vm.path <- system.file("extdata", "MINE.jar", package = "lucky")

  ### 运行MINE函数
  nc.control <- match(control.markers,colnames(expr2))
  nc.target <- setdiff(1:ncol(expr2),nc.control)
  print("批量运行MINE...")
  ## 运行MINE
  library("rJava")
  #i=2
  .jinit(classpath=java.vm.path)
  newcolnames <- c("X var","Y var","MIC (strength)","MIC-p^2 (nonlinearity)","MAS (non-monotonicity)","MEV (functionality)","MCN (complexity)","Pearson(p)")
  if(method == "one.pair"){
    result <- NULL
    for(i in 1:length(nc.target)){
      ## 单对
      MINE("MINE.test.csv",method,nc.control-1,nc.target[i]-1)
      ## 读取Results.csv中的信息
      result.path <- list.files(path = "./",pattern = ",Results.csv",full.names = T);result.path
      result1 <- read.csv(result.path)
      colnames(result1) <- newcolnames
      result <- rbind(result,result1)
      ## 删除多余的文件
      print("正在清除垃圾文件...")
      delete.path <- list.files(path = "./",pattern = "MINE.test.csv,",full.names = T);file.remove(delete.path)
      print("完成垃圾文件清除!")
      }
  } else {
    ## 全部对
    MINE("MINE.test.csv","all.pairs",0)
    ## 读取Results.csv中的信息
    result.path <- list.files(path = "./",pattern = ",Results.csv",full.names = T);result.path
    result <- read.csv(result.path)
    colnames(result) <- newcolnames
    ## 删除多余的文件
    print("正在清除垃圾文件...")
    delete.path <- list.files(path = "./",pattern = "MINE.test.csv,",full.names = T);file.remove(delete.path)
    print("完成垃圾文件清除!")
  }
  ##保存文件
  save(result,file = paste0("result_",method,".rda"))

  ###输出结果
  if(method=="one.pair"){
    return(result)
  } else {
    ### 有所有对的资料，可以输出相关矩阵
    ## result生成配对
    a1 <- as.character(result[,1])
    a2 <- as.character(result[,2])
    names <-unique(a1)

    ### 组装数据为相关矩阵
    mt1 <- matrix(rep(0,length(names)*length(names)),
                  nrow = length(names),
                  dimnames = list(names,names))
    for(i in 1:length(names)){
      for(j in 1:length(names)){
        mt1[i,j] <- paste(rownames(mt1)[i],colnames(mt1)[j],sep = "-")
      }
    } #填上配对
    mt2 <- as.character(mt1)
    #pair.i="Area-Frost"
    getcor <- function(pair.i){
      pair.i1 <- unlist(strsplit(pair.i,"-"))
      l1 <- a1 %in% pair.i1
      l2 <- a2 %in% pair.i1
      l3 <- ifelse(l1==T & l2==T,T,F)
      cor.i <- result$`MIC (strength)`[l3]
      cor.i <- as.numeric(as.character(cor.i))
      return(cor.i[1])
    }
    cor <- apply(as.matrix(mt2),1,getcor)
    print("正在形成相关矩阵...")
    cor1 <- matrix(cor,
                   nrow = length(names),
                   dimnames = list(names,names))
    print("完成相关矩阵组装!")
    cor1[is.na(cor1)] <- 1
    result1 <- list(MINE.result=result,MIC.matirx=cor1)
    return(result1) #View(result1[["MIC.matirx"]])
  }
  #End
}


##Value
# X var/Y var  # 作比较的markers pair.
# MIC (strength) # the maximal information coefficient，相关性强度。MIC越高，相关性越强。与Pearson correlation是类似的。
# MIC-p^2 (nonlinearity) #MIC-Pearson^2，非线性程度。y=x的MIC-p^2为0。MIC-p^2值越高，相关性的非线性程度就越强。
# MAS (non-monotonicity) # 非单调性。y=x的MAS为0。MAS越高，则非单调性越强,越难用单调函数进行相关性描述。
# MEV (functionality) # 函数性。y=x的MEV为1,MEV越高代表函数性越强，越容易用函数进行描述。
# MCN (complexity) #复杂度:y=x的值为2，越高表明相关性越复杂，越难用尽量少的变量去描述它们之间的相关性。
# Pearson(p) #Pearson correlations




