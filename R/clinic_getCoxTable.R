

# getCoxTable help get a model result of univariate and multivariate Cox regression Models.

# data# a data frame containing cluster,time and status value.
# time.col#colnames of time value
# status.col#colnames of status value
# stepwise#Whether use stepwise strategy to simple multiple cox model
# direction="both"#method of stepwise.See MASS::stepAIC.
# cluster#the cluster values like age and sex.
# control # set the control value in a factor vector.For example,in T status,T1 should be a control in the contrast of T3 vs T1.
# dig=2#the decimal place of output numeric value.
# names#part of file name.
#' @export
getCoxTable <- function(data,
                        time.col,
                        status.col,
                        stepwise=F,direction="both",
                        cluster,control,
                        dig=2,
                        names = "test1"){
  ## 包
  nd <- c("survival","broom","MASS","dplyr")
  Plus.library(nd)

  ## 提取数据
  data1 <- data[,c(time.col,status.col,cluster)]
  colnames(data1)[1:2] <- c("time","status")

  ## 重因子化
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

  ### 单因素分析
  Univariate <- NULL;Univariate.ph <- NULL
  for(i in 1:length(cluster)){ # i=1
    c.i <- cluster[i]
    data1.i <- subset(data1,select = c("time","status",c.i))
    fit.i <- coxph(Surv(time,status)~., data = data1.i)
    result.i <- summary(fit.i)[["coefficients"]]
    Univariate <- rbind(Univariate,result.i)
    # ph假定
    ph.i <- cox.zph(fit.i)
    Univariate.ph <- rbind(Univariate.ph,ph.i$table)
  }
  Univariate <- as.data.frame(Univariate)
  Univariate.ph <- as.data.frame(Univariate.ph)

  ### 多因素分析
  ## 是否进行stepwise
  if(stepwise == T){
    ## 除去数据中的空值
    print("运行MASS::stepAIC程序...")
    logic1 <- apply(data1,1,is.one.na)
    data2 <- data1[!logic1,];test1 <- table(logic1)[names(logic1) %in% T]
    if(T %in% logic1){print(paste0("There are ",test1," rows containing at least 1 NA value."))}
    ## cox回归
    fit <- coxph(Surv(time,status) ~.,data = data2)
    ##stepwise
    step1 <- stepAIC(fit,direction=direction)
    fit2 <- coxph(step1[["formula"]],data = data2)
    print("完成逐步AIC筛选法!")
  } else {
    fit2 <- coxph(Surv(time,status) ~.,data = data1)
  }
  m <- summary(fit2)
  multivariate <- m[["coefficients"]]
  multivariate <- as.data.frame(multivariate)
  # multiple ph
  multivariate.ph <- cox.zph(fit2)
  multivariate.ph <- as.data.frame(multivariate.ph$table)

  ### 结果汇总
  # 单因素分析
  ci <- -round(qnorm((1-0.95)/2),2)
  Factor = rownames(Univariate)
  hr= Univariate$`exp(coef)`
  low.hr =  exp(Univariate$coef - ci*Univariate$`se(coef)`)
  up.hr =  exp(Univariate$coef + ci*Univariate$`se(coef)`)
  hr1 = paste(round(hr,dig),"(",round(low.hr,dig)," to ",round(up.hr,dig),")",sep = "")
  p.val = as.character(round2(Univariate$`Pr(>|z|)`,dig))
  df.u <- data.frame(
    Factor=Factor,
    hr=hr1,
    p=p.val
  )

  ## 多因素
  # 单因素分析
  Factor = rownames(multivariate)
  hr= multivariate$`exp(coef)`
  low.hr =  exp(multivariate$coef - ci*multivariate$`se(coef)`)
  up.hr =  exp(multivariate$coef + ci*multivariate$`se(coef)`)
  hr1 = paste(round(hr,dig),"(",round(low.hr,dig)," to ",round(up.hr,dig),")",sep = "")
  p.val = as.character(round2(multivariate$`Pr(>|z|)`,dig))
  df.m <- data.frame(
    Factor=Factor,
    hr=hr1,
    p=p.val
  )

  df1 <- left_join(df.u,df.m,by="Factor")
  df1 <- as.matrix(df1)

  ## 重命名行名
  rn <- df1[,"Factor"]
  newrownames <- function(rn,cluster.i,control.i){
    if(is.na(control.i)){
      #说明名字不需要更改
      name <-  cluster.i
    } else {
      #说明名字要更改
      p1 <- grep(cluster.i,rn)
      rn.i <- rn[p1]
      n2 <- Fastextra(rn.i,cluster.i,2)
      n1 <- rep(cluster.i,length(n2))
      name <- paste(n1,n2,sep = "_")
      name <- paste(name,control.i,sep = " vs. ")
    }
    return(as.character(name))
  }
  test1 <- data.frame(cluster = cluster,control = control)
  new.rn <- apply(test1,1,function(x)newrownames(rn,x[1],x[2]))
  df1[,"Factor"] <- unlist(new.rn)

  ## 进一步处理
  title <- c("Factor","HR(95%CI)","P","HR(95%CI)","P")
  title1 <- c(" ","Univariate analysis"," ","Multivariable analysis"," ")
  df2 <- rbind(title1,title,df1)
  colnames(df2) <- c("Factor","Univariate.HR(95%CI)","Univariate.P","Multivariable.HR(95%CI)","Multivariable.P")

  ##输出结果
  write.csv(df2,paste0(names,"_Cox proportional hazards models.csv"),row.names = F)
  print("Cox模型结果已保存在当前工作目录。")
  l <- list(
    Univariate=Univariate,
    Multivariate = multivariate,
    PHtest = list(Univariate = Univariate.ph,
                  Multivariate =  multivariate.ph)
  )
  return(l)
}





