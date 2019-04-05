

#通常建议用nonrandom的策略来进行匹配。结果输出包括匹配前后的对比、PS分布情况，并最后输出匹配后的design文件。PSM的方法在文献中广泛应用,它可以较好地平衡两组均衡性与样本数。pair.random对于异质性较强的数据难以得到尽可能多的matched pair。但如果对样本量无特殊要求，应该使用pair.random进行精准matched

###========pair.random
#用于对患者按某对contrast进行1:1抽样。
# data = design1
# contrast="N.status"
# cluster = c("age1","gender","T.status")
# seed = 2018
# ratio=0.5
# pair = TRUE
# summary = F
#' @export
pair.random <- function(data,
                        contrast,
                        cluster,
                        seed = 2018,
                        ratio=0.2,
                        pair = TRUE,
                        summary = T){
  library(stringr)
  #记录ID
  data$ID <- rownames(data)

  #筛选数据1：去除NA值。
  vt1 <- c(contrast,cluster)
  p <- NULL
  for (i in vt1) {
    p.i <- match(i,colnames(data))
    p.i <- data[,p.i]
    p.i <- which(is.na(p.i) == T)
    p <- c(p,p.i)
  }
  p <- unique(p)
  if(length(p)>0){data1 <- data[-p,]}else{data1 <- data}
  data1

  ##按新分类来进行区分
  p1 <- match(cluster,colnames(data1))
  data1$col1 <- apply(data1[,p1],1,function(x) paste(x,collapse  = "_"))
  data1 <- data1[order(data1$col1),]

  ##按contrast将data1分开
  l1 <- list();v <- NULL;
  p1 <- match(contrast,colnames(data1))
  level <- as.character(unique(data1[,p1]))
  i = "A"
  for (i in level) {
    p.i <- grep(i,data1[,p1])
    data.i <- data1[p.i,]
    l1 <- c(l1,list(data.i))
  }
  names(l1) <- level

  ##众多新分类中的共同分类
  v <- NULL;
  for(i in 1:length(l1)){
    data2 <- l1[[i]]
    l2.i <- names(table(data2$col1))
    v <- c(v,l2.i)
  }
  can.v <- names(table(v))[table(v)>1]
  cannot.v <- names(table(v))[table(v)==1]
  if(length(length(can.v)==length(v))){
    print(str_c("seed_",seed,":所有分类都可进行分别抽样。"))
  } else {
    print(paste("seed_",seed,":由于并非共同分类，以下分类无法在contrast中分别抽样：",paste(cannot.v,collapse = ";")))
  }

  if(pair == T){
    #按分类含量较少的标准来进行两组的选择。
    v2 <- NULL
    for(i in 1:length(l1)){
      data2 <- l1[[i]]
      data3 <- NULL;
      for(j in can.v){
        p.j <- grep(j,data2$col1)
        data3.j <- data2[p.j,]
        data3 <- rbind(data3,data3.j)
      }
      v2.i <- table(data3$col1)
      v2 <- c(v2,v2.i)
    }
    min1 <- NULL;
    for(i in unique(names(v2))){
      min1.i <- min(v2[names(v2)==i])
      min1 <- c(min1,min1.i)
    }
    names(min1) <-unique(names(v2))
    #提取某分类的id
    extra.id <- function(data,sc,seed,ratio){
      select <- NULL;
      for(i in sc){
        data.i <- data[data$col1 == i,]
        c.p <- match(i,names(min1))
        c <- min1[c.p]
        set.seed(seed);select.i <- sample(data.i$ID,floor(ratio*c))
        select <- c(select,select.i)
      }
      return(select)
    }
    v1 = NULL;
    for(i in 1:length(l1)){
      v1.i <- extra.id(l1[[i]],can.v,seed,ratio)
      v1 <- c(v1,v1.i)
    }
    print("random pair completed")
  } else {
    #不按pair来取值
    l2 <- list()
    for(i in 1:length(l1)){
      data2 <- l1[[i]]
      data3 <- NULL;
      for(j in can.v){
        p.j <- grep(j,data2$col1)
        data3.j <- data2[p.j,]
        data3 <- rbind(data3,data3.j)
      }
      l2.i <- table(data3$col1)
      l2<- c(l2,list(l2.i))
    }
    #提取某分类的id
    extra.id <- function(data,l2i,seed,ratio){
      select <- NULL;
      for(i in names(l2i)){
        data.i <- data[data$col1 == i,]
        c.p <- match(i,names(l2i))
        c <- l2i[c.p]
        set.seed(seed);select.i <- sample(data.i$ID,floor(ratio*c))
        select <- c(select,select.i)
      }
      return(select)
    }
    v1 = NULL;
    for(i in 1:length(l1)){
      v1.i <- extra.id(l1[[i]],l2[[i]],seed,ratio)
      v1 <- c(v1,v1.i)
    }
    print("random pair completed")
  }

  ##order
  x1 <- data[match(v1,data$ID),]
  x2 <- as.matrix(x1)
  order.col <- NULL;
  for(a in cluster){
    order.col.i <- grep(a,colnames(x2))
    order.col <- c(order.col,order.col.i)
  }
  query <- paste('x2[,',order.col,']',sep='') %>%
    paste(.,collapse=',') %>%
    paste('order(',.,',decreasing =T)',sep='')
  class(query);x2 <- x2[eval(parse(text = query)),]
  x2 <- as.data.frame(x2)

  ##选择后的data
  if(summary == T){
    x2 <- subset(x2,select = c(cluster,contrast))
    return(x2)
  } else {
    library(plyr)
    p <-  match("ID",colnames(x1))
    x2 <- x2[,-p]
    return(x2)
  }
  #End
}




