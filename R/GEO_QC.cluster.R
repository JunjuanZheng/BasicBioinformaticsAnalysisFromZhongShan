


###======QC.cluster()
# QC.cluster()通过选取少数基因进行聚类，从而判断分组的差异程度。QC.cluster()基于hwb_BI::heatmap.dds()函数。输出结果为heatmap图。
#' @export
QC.cluster <- function(DesignEset,
                       contrast=NULL,
                       select.c=1000,#选取的基因数
                       seed=2018,
                       savefile=F,
                       log.convert = F,
                       names="GSE1231",#part of file name
                       cluster_rows = T,
                       cluster_cols =F,
                       row.k = 4){
  ##加载包
  library(stringr);library(limma);library(Biobase)

  ## contrast的设置
  if(is.null(contrast)){
    #如果不提供contrast，就直接用全部的contrast
    contrast.matrix <- DesignEset[["DesignList"]][["contrast.matrix"]]
    c1 <- colnames(contrast.matrix)
    c2 <- list()
    for(i in c1){
      c2.i <- unlist(strsplit(i,"-"))
      c2 <- c(c2,list(c2.i))
    }
  } else {
    #提供了contrast，直接使用定义的contrast
    c2 <- contrast
  }

  ##提取某个DesignList和contrast
  p1 <- NULL #i=1
  for(i in 1:length(c2)){
    ## 提取某个DesignList和contrast
    contrast.i <- c2[[i]]
    DesignEset.i <- getOneDesignEset(DesignEset,contrast.i)
    names.i <- paste0(contrast.i,collapse = "-")

    ## condition
    cond <- DesignEset.i[["DesignList"]][["condition"]]

    ## 表达矩阵
    library(limma);library(Biobase)
    expr <- exprs(DesignEset.i[["ESet"]])

    ## 将control探针排除在外(从QC.cluster内移植过来)
    annotation0 <- fData(DesignEset.i[["ESet"]])
    control.probes <- DesignEset.i[["DesignList"]][["control.probes"]]
    if(is.null(control.probes)){
      #并无control.probes
      expr1 <- expr
    } else {
      #有对照探针
      con.position <- NULL
      for(i in 1:length(control.probes)){
        con.i <- control.probes[[i]]
        #找到control.probes的位置
        con.position.i <- NULL
        # s="control"
        for(s in con.i$control.symbol){
          con.position.s <- grep(s,annotation0[,con.i$control.col])
          con.position.i <- c(con.position.i,con.position.s)
        }
        con.position <- c(con.position,con.position.i)
      }
      con.position <- unique(con.position)
      # 去除control探针
      all.position <- 1:length(annotation0[,con.i$control.col])
      l1 <- all.position %in% con.position
      expr1 <- expr[!l1,]
    }
    select1 <- rownames(expr1)

    # select.c:选取部分基因
    set.seed(seed);select.c1 <- sample(1:nrow(expr1),select.c)
    if(is.null(select.c1)){
      select2=select1
    } else {
      select2=select1[select.c1]
    }

    ## heatmap
    # contrast.list
    set3.order <- c(1,4,6,10,2,3,5,7,8,9,11,12)
    library(RColorBrewer)
    c <- c("red","green",brewer.pal(12,"Set3")[set3.order])
    contrast.list = list(group=contrast.i,color=c[1:length(contrast.i)])
    p.i <- heatmap.dds(expr.matrix =  expr1,
                      select = select2,
                      design = cond,
                      contrast = "condition",
                      contrast.list=contrast.list,
                      rowscale=T,log.convert =log.convert,
                      expr.name = "expression",
                      cluster_rows = cluster_rows,
                      cluster_cols=cluster_cols,
                      row.k = row.k,
                      show_column_names = F,
                      show_row_names = F)
    p1 <- c(p1,list(p.i))

    ## 预览图片
    win.graph(width = 9,height = 8);print(p.i)
  }

  #保存文件
  if(savefile==T){
    pdf(str_c(names,"_QualityHeatmap.pdf"),width=9,height=8)
    for(i in 1:length(p1)){
      p.i <- p1[[i]]
      print(p.i)
    }
    dev.off()
    print(str_c(names,"_QualityHeatmap.pdf已保存在本地"))
    #End
    return(p1)
  } else {
    #End
    return(p1)
  }

}

# data(DesignEset)
# contrast=list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")) # a list that contain lots of contrasts.
# QC.cluster(DesignEset,
#            contrast,
#            select.c=1000,#选取的基因数
#            seed=2018,
#            savefile=F,
#            log.convert = F,
#            names="GSE1231",#part of file name
#            cluster_rows = T,
#            cluster_cols =F,
#            row.k = 4)


