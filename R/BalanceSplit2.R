



BalanceSplit2 <- function(data,
                          y,
                          cluster,
                          nsplit = 2,
                          seed.range=1:500,
                          p.cutoff=0.1){
  ## 加载必要的包
  need <- c("caret","knitr","tableone")
  Plus.library(need)

  ## 给出y值的数据，按y值排序
  y = data[,y]

  ## 生成模型列表
  L1 <- NULL;
  for(z in 1:length(seed.range)){ #z=1
    seed.i <- seed.range[z]
    set.seed(seed.i);
    NO.list <- createDataPartition(y = y, times = nsplit)


    ## 构建数据集
    data.list <- NO.list
    for(i in 1:length(NO.list1)){
      select.i <- NO.list1[[i]]

      data.i <- data[select.i,]
      data.list[[i]] <- data.i
      names(data.list)[i] <- i
    }

    ## 组装数据集
    data2 <- NULL
    for(i in 1:length(data.list)){
      data.i <- data.list[[i]]
      data.i$Set <- as.character(i)
      data2 <- rbind(data2,data.i)
    }

    ## 获得某个数据集的行名
    cohort <- NULL
    for(i in 1:length(data.list)){
      cohort.i <- rownames(data.list[[i]])
      cohort <- c(cohort,list(cohort.i))
      names(cohort)[i] <- names(data.list)[i]
    }

    ## TableOne
    tableOne <- CreateTableOne(vars = cluster, strata = "Set", data = data2)
    table1 <- print(tableOne,
                    printToggle = FALSE,
                    noSpaces = TRUE)

    ## P值处理
    p.val <- as.character(table1[,"p"])
    p.val <- as.numeric(as.character(p.val))
    p.val <- p.val[!is.na(p.val)]

    ## 筛选
    if(all(p.val>p.cutoff)){
      #p值全部大于0.05
      print(paste0("seed",seed.i,"满足条件!"))
      l.i <- list(seed=seed.i,
                  p.min=min(p.val),
                  NO.list=NO.list1,
                  data=data2,
                  table=table1,
                  cohort = cohort)
      L1 <- c(L1,list(l.i))
      names(L1)[length(L1)] <- paste0("seed_",seed.i)
    } else {
      median=NULL
    }
  }







}












