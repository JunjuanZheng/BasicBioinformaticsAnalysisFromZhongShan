

####======================get.journal.ISSN=======================####
#按杂志名称自动获取ISSN号（网络版）。杂志较多是比较耗时。程序结束可选择自动播放音乐。
get.journal.ISSN <- function(journal.names,show.music = T){

  #加载包
  library(pacman)
  p_load(rvest);p_load(tuneR)
  library(stringr);library(stringi)

  #report1
  print(str_c("共有",length(journal.names),"本杂志的ISSN号待搜索..."))

  #按URL获得ISSN号
  df <- NULL;
  for(n in 1:length(journal.names)){
    #转换格式
    #journal.i <- journal.names[n]
    journal.i <- Fastextra(journal.names[n],"-",1)
    special1 = c("[ ]","-","&","_");
    special2 = c("+","-","%26","_")
    for(i in 1:length(special1)){
      journal.i <- gsub(special1[i], special2[i], journal.i)
    }

    gurl <- str_c("https://portal.issn.org/api/search?search[]=MUST=keyproper,keyqualinf,keytitle,notcanc,notinc,notissn,notissnl,unirsrc=",journal.i)
    a <- tryCatch(md <- gurl %>%  read_html(encoding="GBK"), error = function(e)e, finally = NULL)
    if(is.null(a[["message"]])){
      ##程序没有报错
      #提取ISSN号
      ISSN <- md %>%
        html_nodes("div.squaredFour") %>%
        html_nodes("input") %>%
        html_attr("value")

      if(length(ISSN) == 0){
        print(str_c("无法正确提取第",n,"个杂志的ISSN编码。"))
        #说明url虽然没错，但是无法获取ISSN
        df.i <- data.frame(journal =journal.names[n],search.names = "ISSN.NULL",ISSN.code ="ISSN.NULL")
      } else {
        #提取名称
        jn <- md %>%
          html_nodes("div.item-result-block-title") %>%
          html_nodes("h5") %>%
          html_text()

        #提取国家
        #c <- md %>%
        # html_nodes("div.flex-zero") %>%
        # html_nodes("p") %>%
        #html_text()
        #c1 <- c[c(2,7,11)]
        #c2 <- unlist(strsplit(c1,"[:]"))[c(2,4,6)]

        ##组合
        p <- grep("Online",jn)
        if(length(p) == 0){
          print(str_c("第",n,"个杂志并没有网络版"))
          #说明并无Online的版本。不适合用于网络检索。
          df.i <- data.frame(journal =journal.names[n],search.names = "NON.Online",ISSN.code =paste0(ISSN,collapse = "_"))
        } else {
          #说明有Online的版本，适合用于网络检索。
          df.i <- data.frame(journal =journal.names[n],search.names = jn[p],ISSN.code =ISSN[p])
          print(str_c("完成第",n,"个杂志的ISSN编码搜索。"))
        }
      }

    } else {
      ##程序有报错。
      print(str_c("地址解析错误，无法完成第",n,"个杂志的ISSN编码搜索。"))
      ##组合
      df.i <- data.frame(journal =journal.names[n],search.names = "error",ISSN.code ="error")
    }
    df <- rbind(df,df.i)
  }

  #输出结果
  if(show.music == T){
    tuneR::play(music)
    return(df)
  } else {
    return(df)
  }

}


####============================get.ISSN=========================##
#'after get.journal.ISSN
#'@description get.ISSN is after get.journal.ISSN and produce unique ISSN code via intelligent algorithm.
#'@param data a data containing journal,ISSN.code and search.names,which is always produce by luckey::get.journal.ISSN
#'@param journal.name.col the colname of journals
#'@param ISSN.col the colname of ISSN codes
#'@param search.col the colname of search names
#'@author Weibin Huang<\email{654751191@@qq.com}>
#'@examples
#'## get unique ISSN
#'library(lucky)
#'df1 <- get.ISSN(ISSN.3.Inf)
#'df2 <- get.ISSN2(df1)
#'
#'## get IF information
#'data("IF")
#'colnames(df2)[1] <- "Full Journal Title"
#'library(dplyr)
#'ISSN_3.Inf <- left_join(df2,IF,by = "Full Journal Title")
#'table(is.na(ISSN_3.Inf$`Journal Impact Factor`))
#'
#'## get new IF plus ISSN
#'colnames(ISSN_3.Inf) <- c("Journals","Search names","ISSN","Total #'Cites","IF","Eigenfactor Score")
#'ISSN_3.Inf <- subset(ISSN_3.Inf,select = c("Journals","ISSN","IF","Search names","Total Cites","Eigenfactor Score"))
#'IF_3.Inf <- ISSN_3.Inf
#'save(IF_3.Inf,file = "IF_3.Inf.rda")
#'@export
get.ISSN <- function(data,
                     journal.name.col="journal",
                     ISSN.col="ISSN.code",
                     search.col="search.names"){

  ## 加载必要的包
  need <- c("plyr","stringi")
  Plus.library(need)

  ## 提取journal names
  jn <- factor(unique(data[,journal.name.col]))

  ## 对于每个杂志进行特别处理 jn.i = jn[11]
  ISSN.1 <- function(jn.i){
    # 获得杂志名
    jn.i <- as.character(jn.i)

    # 提取数据
    df1 <- data[data[,journal.name.col] %in% jn.i,]

    # 处理数据
    if(nrow(df1)==1){
      #记录唯一，直接输出结果
      df2 <- df1
    } else {
      #记录不唯一，要选择正确结果
      a <- tolower(jn.i)
      A <- tolower(as.character(df1[,search.col]))
      A <- Fastextra(A," [(]",1)
      logic1 <- A %in% a
      if(T %in% logic1){
        p.i <- match(T,logic1)
        df2 <- df1[p.i,]
      } else {
        df2 <- df1[1,]
      }
    }

    # 输出结果
    return(df2)
  }

  ## 批量操作
  print("an intelligent process is running,please wait...")
  df1 <- adply(as.matrix(jn),1,ISSN.1)
  print("intelligent process done!")
  #table(is.na(df1$journal))

  ## 输出结果
  df2 <- df1[,-1]
  return(df2)

}
















