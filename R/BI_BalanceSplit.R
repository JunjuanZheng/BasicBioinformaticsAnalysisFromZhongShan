
# BalanceSplit help quickly splitting a data frame into balanced parts according to specified cluster factor.

# data # a data frame containing cluster cols.
# cluster# cluster colnames.
# nsplit# the number of balance splitting.Default is 2.
# seed.range# candidate seeds.Default is 1:500.
# p.cutoff# the cut off of P value.Default is 0.1.
#' @export
BalanceSplit <- function(data,
                         cluster,
                         nsplit = 2,
                         seed.range=1:500,
                         p.cutoff=0.1){
  ## 加载必要的包
  need <- c("knitr","tableone")
  Plus.library(need)

  ## 生成模型列表
  L1 <- NULL;
  for(z in 1:length(seed.range)){ #z=1
    seed.i <- seed.range[z]

    ## 随机种子产生
    c <- 1:nrow(data)
    NO.list <- cut.vector(c,nsplit = nsplit)
    NO.list1 <- NO.list
    for(i in 1:length(NO.list)){
      c.i <- length(NO.list[[i]])
      set.seed(seed.i);c.i2 <- sample(c,c.i);c.i2
      NO.list1[[i]] <- c.i2
      c <- setdiff(c,c.i2)
    }
    #intersect(NO.list1[[1]],NO.list1[[2]])#无交集

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

  ## 排序
  print(paste0("完成筛选!共有",length(L1),"个种子符合条件。正在排序..."))
  p.min <- NULL
  for(z in 1:length(L1)){
    L.i <- L1[[z]]
    p.min <- c(p.min,L.i[["p.min"]])
  }
  names(p.min) <- names(L1)
  p.min1 <- sort(p.min,decreasing = T)
  p1 <- Fastmatch(names(p.min1),names(p.min))
  L2 <- L1[p1]

  ## 输出结果
  print("完成排序!")
  return(L2)
}


summary.BalanceSplit <- function(BalanceList,
                                 upset.seed=c(3,11),
                                 set=1){
  ## upset.n
  L2 <- BalanceList
  upset.seed1 <- paste("seed_",upset.seed,sep = "")
  upset.n <- Fastmatch(upset.seed1,names(L2))

  ## get list
  L3 <- L2[upset.n]
  ul <- NULL
  nsplit <- length(unique(L2[[1]]$data$Set))
  for( c in 1:nsplit){# c=1
    ul.c <- NULL
    for(i in 1:length(L3)){#i = 1
      ul.c.i <- L3[[i]]$cohort[[c]]
      ul.c <- c(ul.c,list(ul.c.i))
      names(ul.c)[i] <- names(L3)[i]
    }
    ul <- c(ul,list(ul.c))
    names(ul)[c] <- c
  }

  ## select a set
  ul.c <- ul[[set]]
  win.graph(width = 12,height = 8);
  a <- FastUpset(data = ul.c,
            mainbar.y.label = "Intersection Size",
            sets.x.label = "Set Size",
            order.by = c("freq","degree")[2])
  return(a)
}





