

# Optimized baselined information production from tableone package::CreateTableOne
# CreateTableOne2 is a simplified and enhanced version of table::CreateTableOne.It can refactor the data frame and make visualization and file saving easier.

# cluster# the row markers like Gender
# control # the control of row marker like FEMALE
# strata # the col marker like N.status
# data# a data frame
# save.file # whether to save file
# names# part name of saved file
#' @export
CreateTableOne2 <- function(cluster,
                            control = NULL,
                            strata,
                            data,
                            save.file = T,
                            names = "love"){
  ## 加载必要的包
  nd <- c("tableone","knitr")
  Plus.library(nd)

  ## 选择数据
  data1 <- subset(data,select = c(strata,cluster))

  ## 重因子化
  if(!is.null(control)){
    #对因子顺序有要求
    for(i in 1:length(cluster)){
      control.i <- control[i]
      u.i <- unique(as.character(data1[,cluster[i]]))
      if(is.na(control.i)){
        #连续型变量，不需要因子化
        data1[,cluster[i]] <- data1[,cluster[i]]
      } else {
        #重因子化
        data1[,cluster[i]] <- factor(data1[,cluster[i]],levels = c(control.i,setdiff(u.i,control.i)))
      }
    }
  } else {
    #对因子顺序无要求
    for(i in 1:ncol(data1)){
      if(!is.null(levels(data1[,i]))){
        #有levels的factor变量
        data1[,i] <- factor(data1[,i],levels = as.character(unique(data1[,i])))
      }
    }
  }

  ## tableone
  table1 <- CreateTableOne(vars = cluster,
                           strata = strata,
                           data = data1)
  table1 <- print(table1,printToggle = FALSE,noSpaces = TRUE)
  n <- length(as.character(unique(strata))) + 1
  t1 <- kable(table1[,1:n],align = 'c')
  print(t1)

  ## 保存文件
  if(save.file == T){
    write.csv(table1,paste0(names,"_baseline information.csv"))
  }

  ## 输出结果
  return(t1)
}







