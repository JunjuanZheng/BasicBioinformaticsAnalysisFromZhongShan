

###======pair.PSM
#通过倾向性评分对样本进行1：1的matched。支持"nonrandom","MatchIt"两种策略。
#MatchIt::matchit.method = "exact" (exact matching), "full" (full matching), "genetic" (genetic matching), "nearest" (nearest neighbor matching), "optimal" (optimal matching), and "subclass" (subclassification) are available. The default is "nearest"
#' @export
pair.PSM <- function(data,
                     contrast,
                     treat,
                     cluster ,
                     strategy,
                     method,
                     seed = 2018,
                     ratio=1,
                     summary = F){
  ##下载或加载包
  library(pacman)
  p_load(wakefield);p_load(tableone);p_load(knitr);p_load(stringr)

  ##记录id
  data.raw <- data
  data <- as.data.frame(data)
  data$ID <- rownames(data)

  ##重新定义因子型变量的level，不然在进行OneTable运算时容易有bug。
  for(i in cluster){
    if(is.factor(data[,i])){
     data[,i] <- factor(data[,i],levels = as.character(unique(data[,i])))
    } else {
      print(str_c(i," is not a factor value."))
    }
  }

  ##去除变量中的NA值
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
  rownames(data1) <- 1:nrow(data1)

  ##初步展示数据匹配前的异质性
  table1 <- CreateTableOne(vars = cluster,
                           data = data1,
                           strata = contrast)
  table1 <- print(table1,
                  printToggle = FALSE,
                  noSpaces = TRUE)
  b <- kable(table1[,1:3],align = 'c')

  ##PSM计算：有两种方法可以选，分别来自MatchIt包和nonrandom包
  if(strategy == "MatchIt"){
    ###以MatchIt的方法来进行PSM的计算
    p_load(MatchIt);p_load(optmatch)
    #logistic回归公式构建
    data2 <- subset(data1,select = c("ID",contrast,cluster))
    data2$Group <- ifelse(data2$N.status == treat,TRUE,FALSE)
    f <- paste0("Group ~ ",paste0(cluster,collapse = " + "))
    f <- as.formula(f)
    #PSM计算
    set.seed(seed)
    match.it <- matchit(formula = f,
                        data = data2,
                        method = method,#这个要按要求变化
                        ratio = ratio)
    match.a <- summary(match.it)
    a <- kable(match.a$nn,
               digits = 2,
               align = 'c',
               caption = 'Table 2: Sample sizes')
    #匹配后的年龄和性别分布基本一致了
    a2 <- kable(match.a$sum.matched[c(1,2,4)],
                digits = 2,
                align = 'c',
                caption = 'Table 3: Summary of balance for matched data')

    #结果输出
    pair <- match.data(match.it)[1:ncol(data2)]
    print("匹配前的情况：");print(b)
    print("匹配后的情况：");print(a2)
    print(plot(match.it,type = 'jitter',interactive = FALSE))
  } else {
    if(strategy == "nonrandom"){
      ###以nonrandom包来进行PSM
      p_load(nonrandom)
      data2 <- subset(data1,select = c("ID",contrast,cluster))
      data2$group = ifelse(data2[,contrast] == treat,1,0)
      f <- paste0("group ~ ",paste0(cluster,collapse = " + "))
      f <- as.formula(f)
      #PSM计算
      match.it2 <- pscore(data = data2,formula = f)
      par(mfrow = c(1,2))
      print("匹配前：")
      print(b)
      plot.pscore(x = match.it2,
                  main = "Pre-matched PS distribution",
                  xlab = "",
                  par.1=list(col="red"),#处理组是红色
                  par.0=list(lwd=2),
                  par.dens=list(kernel="gaussian"))
      logic1 <- ifelse(length(which(data2$group==1)) > length(which(data2$group==0)),F,T)#treat匹配至untreated。一般treat少的时候，应该设置为FALSE
      list.match <- ps.match(object = match.it2,
                             who.treated =1,#表示1代表处理组
                             ratio = ratio,#匹配比例：1：ratio
                             caliper = "logit",
                             x = 0.2,#得分容差?
                             givenTmatchingC = logic1,
                             matched.by= "pscore",
                             setseed = seed)
      #summary(list.match)
      pair <- list.match$data.matched

      #查看PSM的效果
      stable1 <- CreateTableOne(vars=cluster,
                                strata=contrast,
                                data=pair)
      #print(stable1,showAllLevels = TRUE)
      stable1 <- print(stable1,
                       printToggle = F,
                       noSpaces = TRUE)
      a <- kable(stable1[,1:3],align = 'c')
      print("匹配后：")
      print(a)
      match.it3 <- pscore(data = pair[,-match("pscore",colnames(pair))],
                          formula = f)
      plot.pscore(x = match.it3,
                  main = "After-matched PS distribution",
                  xlab = "",
                  par.1=list(col="red"),#处理组是红色
                  par.0=list(lwd=2),
                  par.dens=list(kernel="gaussian"))
    } else {print("strategy错误，请检查")}

  }

  #End
  par(mfrow = c(1,1))

  ##输出结果
  if(summary == T){data2 <- pair} else {data2 <- data.raw[as.character(pair$ID),]}

  ##order
  #x2 <- as.matrix(data2)
  #order.col <- NULL;
  #for(a in cluster){
  #  order.col.i <- grep(a,colnames(x2))
  #  order.col <- c(order.col,order.col.i)
  #}
  #p_load(dplyr)
  #query <- paste('x2[,',order.col,']',sep='') %>%
  #  paste(.,collapse=',') %>%
  #  paste('order(',.,',decreasing =T)',sep='')
  #class(query);x2 <- x2[eval(parse(text = query)),]
  #data2 <- as.data.frame(x2)

  #返回值
  return(data2)

}





