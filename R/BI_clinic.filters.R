

####=====================clinic.filters=====================####
# 用于按某些条件对临床数据进行过滤。
# data：待分类抽样的数据，
# filter：数据的过滤条件
# filter.parameter：过滤条件的参数。这里支持数值型/因子型过滤条件的<、>、<=、>=、=五种运算。目前的<=、>=运算似乎有bug。待修复。
#' @export
clinic.filters <- function(data,
                           filter,
                           filter.parameter,
                           summary = T){
  #筛选数据1：去除NA值。
  p <- NULL
  for (i in filter) {
    p.i <- match(i,colnames(data))
    p.i <- data[,p.i]
    p.i <- which(is.na(p.i) == T)
    p <- c(p,p.i)
  }
  p <- unique(p)
  if(length(p)>0){data1 <- data[-p,]}else{data1 <- data}


  #筛选数据2：按filter及参数去除不合条件者。
  formula1 <- function(data1,filter,filter.parameter){
    operation1 <- c("<","<=",">",">=","=")
    operation.matrix <- matrix(rep(0,10),ncol = 2)
    extra1 <- function(vt1,split){
      a <- unlist(strsplit(as.character(vt1),split))[2]
      if(is.na(as.numeric(a)) == T) return(a) else return(as.numeric(a))
    }
    #找到某filter的向量
    vt.filter <- data1[,filter]
    w = 1;while (w <= 5){
      operation.matrix[w,] <- c(operation1[w],extra1(filter.parameter,operation1[w]))
      w = w + 1
    }
    p <- which(is.na(operation.matrix[,2])==F)#找到对应的运算符
    #构建过滤逻辑向量
    operation2 <- function(vt,s,o){
      if(o == "<") {return(vt<s)} else {if(o == "<=") {return(vt<=s)} else {if(o == ">") {return(vt>s)} else {if(o == ">=") {return(vt>=s)} else {if(o == "=") {return(vt==s)}}}}}}
    s = ifelse(is.na(as.numeric(operation.matrix[p,2])) == T,operation.matrix[p,2],as.numeric(operation.matrix[p,2]))
    o = operation.matrix[p,1]
    return(operation2(vt.filter,s,o))
  }
  data.filter2 <- data1;data.filter2[,1:ncol(data.filter2)] <- NULL;
  for (i in 1:length(filter)) {
    data.filter2[,i] <- formula1(data1,filter[i],filter.parameter[i])
  }
  data.filter21 <- apply(data.filter2,1,function(x) ifelse(FALSE %in% x,FALSE,TRUE))
  data.filter2$logic <- data.filter21
  data2 <- data1[data.filter21,]

  #提取数据
  if(summary == T){
    data3 <- subset(data2,select = filter)
    return(data3)
  } else {
    return(data2)
  }
}

#Test
# data(data.selfprocess2)
# colnames(data.selfprocess2)
# data.selfprocess2$LN.sum <- as.numeric(as.character(data.selfprocess2$LN.sum))
# data.selfprocess2$OS.time <- as.numeric(as.character(data.selfprocess2$OS.time))
# cf1 <- clinic.filters(data = data.selfprocess2,
#                       filter = c("LN.sum","OS.time"),
#                       filter.parameter = c(">15",">179"),
#                       summary = F)
# cf2 <- clinic.filters(data = data.selfprocess2,
#                       filter = c("LN.sum","OS.time","BMI"),
#                       filter.parameter = c(">15",">179",">15"),
#                       summary = F)
# View(cf1);View(cf2)

