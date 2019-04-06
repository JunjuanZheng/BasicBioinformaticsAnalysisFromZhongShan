

####=====================Annotate.Journal=====================####
#' @title Annotate journal id via crawler
#' @description AnnotateJournal help annotate journal id via letpub website.AnnotateJournal is more powerful than get.ISSN series and is recommanded to get information of journals like ISSN or IF.
#' @keywords AnnotateJournal
#' @param journal.id a vector of journal id
#' @param source the data source.Default is letpub.
#' @param report whether report working information
#' @details The meaning of results from \code{AnnoteJournal} is as following: \cr journal.id: the ids of journals in letpub\cr journals:the names of journals\cr abbr: the abbreviation of journal names \cr issn:the ISSN code of journals \cr IF: the impact factor of journals \cr IF5: the merge impact factor of journals in recent 5 years \cr selfcited: self-cited ratio \cr website: the website of the journal \cr OA: whether a OA journal \cr research: the research direction \cr area: the area of the journal like "USA" \cr regular: the period of article checking \cr BasicSubject: the basic subject of the journal \cr SecondSubject: the second subject of the journal \cr TopJournal: whether a top journal \cr ReviewJournal: whether a review journal \cr speed: the speed of article checking \cr recieve: the ratio of article recieved
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @seealso \code{\link{MakeENlist}}
#' @examples
#' AnnotateJournal(journal.id=1:13000,
#'                 source = "letpub")
#' @export
AnnotateJournal <- function(journal.id=1:13000,
                           source = "letpub",
                           report = T){

  ##加载必要的包
  need <- c("rvest","stringr","plyr")
  Plus.library(need)

  ##提取一个id的函数
  ano1 <- function(journal.id,
                   source = "letpub",
                   report = T){
    # journal.id = 2

    ## source
    if(source == "letpub"){
      #构建url
      gurl <- str_c("http://www.letpub.com.cn/index.php?journalid=",journal.id,"&page=journalapp&view=detail")

      ## 错误控制
      a <- tryCatch(md <- gurl %>%  read_html(encoding="UTF-8"), error = function(e)e, finally = NULL)

      if(is.null(a[["message"]])){
        text <- md %>% html_nodes("td") %>% html_text()

        ## 获取杂志名
        title <- md %>%
          html_nodes("td") %>%
          html_nodes("span") %>%
          html_nodes("a")%>%
          html_text();
        title <- title[1];title

        ## 缩写
        subtitle <-md %>%
          html_nodes("td") %>%
          html_nodes("span") %>%
          html_nodes("font") %>%
          html_text();
        subtitle <- subtitle[1];subtitle


        ## ISSN号
        issn <- text
        logic1 <- str_extract(issn,"[0-9]{4}[-][0-9|A-Z]{4}\\b")
        logic2 <- str_extract(issn,"^[0-9]{4}[-][0-9|A-Z]{4}")
        logic <- apply(cbind(logic1,logic2),1,is.one.na)
        issn <- issn[logic==F];issn

        ## IF
        IF <- text
        p1 <- grep("最新IF",IF)[1] #p1 <- grep("最新影响因子",IF)
        p2 <- str_detect(IF,"[0-9][.][0-9]{3}\\b|[0-9]{2}[.][0-9]{3}\\b|[0-9]{3}[.][0-9]{3}\\b")
        p2 <- grep(T,p2)
        if((p1+1) %in% p2){
          #成功找到了IF
          a <- str_extract(IF,"[0-9][.][0-9]{3}\\b|[0-9]{2}[.][0-9]{3}\\b|[0-9]{3}[.][0-9]{3}\\b")
          IF<- a[p1+1]
        } else {
          IF <- "Not Available"
        };IF

        #五年影响因子
        IF5 <- text
        p1 <- grep("五年IF",IF5)+1 #p1 <- grep("五年影响因子",IF5)+1
        IF5 <- IF5[p1];IF5

        ## 自引率
        selfuse <- text
        p1 <- grep("自引率",selfuse)
        selfuse <- selfuse[p1+1]
        selfuse <- selfuse[grep("%",selfuse)[1]];selfuse

        ## h-index
        h_index <- text
        p1 <- grep("h-index",h_index)
        h_index <- h_index[p1+1];h_index

        ## cite score
        cite_score <- text
        p1 <- grep("CiteScore", cite_score)
        cite.score <- cite_score[p1[1]+1]
        cite.score1 <- Fastextra(cite.score,"CiteScore排名",2)
        cite.score2 <- Fastextra(cite.score1,"学科排名百分位")
        score <- str_extract(cite.score2[1],"[0-9]{1,6}[.][0-9]{1,6}")
        q4 <- Fastextra(cite.score2[1],score,2)
        rank <- str_extract(cite.score2[2],"[0-9]{1,6}[ ][/][ ][0-9]{1,6}")
        superior <- Fastextra(Fastextra(cite.score2[2],"小类：",1),"大类：")[2]
        inferior <- Fastextra(Fastextra(cite.score2[2],rank,1),"小类：")[2]


        ## 期刊官方网站
        http <- text
        p1 <- grep("期刊官方网站",http)
        http <- http[ p1+1];http

        ## 是否OA开放访问
        oa <- text
        p1 <- grep("是否OA开放访问",oa)
        oa <- oa[p1+1];oa

        ## 研究方向
        r <- text
        p1 <- grep("涉及的研究方向",r)
        r <- r[p1+1];r

        ## 出版国家或地区
        area <- text
        p1 <- grep("出版国家或地区",area)
        area <- area[p1+1];area

        ##出版周期
        regular <- text
        p1 <- grep("出版周期",regular)
        regular <- regular[p1+1][1];regular

        ## 中科院SCI期刊分区:"大类学科" "小类学科" "Top期刊"  "综述期刊"
        category.status <- md %>%
          html_nodes("td") %>%
          html_nodes("tr") %>%
          html_text()
        category.status <- category.status[grep("区",category.status)][1]
        a <- Fastextra(category.status,"区",3)
        a <- Fastextra(a,"")
        b <- Fastextra(category.status,"区",1:2)
        b <- paste(b,"区",sep="")
        category.status <- c(b,a);category.status

        ##平均审稿速度
        speed <- text
        p1 <- grep("平均审稿速度",speed)
        speed <- speed[p1+1]
        speed <- Fastextra(speed,"：",2);speed

        ## 平均录用比例
        recieve <- text
        p1 <- grep("平均录用比例",recieve)
        recieve <- recieve[p1+1]
        if(length(recieve)==0){
          recieve <- ""
        } else {
          recieve <- Fastextra(recieve,"经验：",2);
          recieve1 <- Fastextra(recieve,"来源Elsevier官网：",1)
          recieve2 <- Fastextra(recieve,"来源Elsevier官网：",2)
          recieve <- paste(recieve1,recieve2,sep = ";");
          if(is.na(recieve2)){recieve <- recieve1}
        }
        recieve

        ## 防止无值产生
        real <- function(x) {
          if(length(x)==0){
            x <- NA
          } else {
            x <- x
          }
          return(x)
        }


        ## 组合成数据框
        if(length(title)==0|is.na(title)){
          #title不存在，即没有对应编号的杂志
          LuckyVerbose(paste0("Error State:Not Such ID:",journal.id,"!"))
          df1 <- data.frame(
            journal.id=journal.id,
            journals="NOT SUCH ID",
            abbr = "NOT SUCH ID",
            issn = "NOT SUCH ID",
            IF = "NOT SUCH ID",
            IF5 = "NOT SUCH ID",
            selfcited = "NOT SUCH ID",
            website = "NOT SUCH ID",
            OA = "NOT SUCH ID",
            research = "NOT SUCH ID",
            area  = "NOT SUCH ID" ,
            regular = "NOT SUCH ID",
            BasicSubject = "NOT SUCH ID",
            SecondSubject = "NOT SUCH ID",
            TopJournal = "NOT SUCH ID",
            ReviewJournal = "NOT SUCH ID",
            speed = "NOT SUCH ID",
            recieve = "NOT SUCH ID"
          )
        } else {
          if(report==T){
            LuckyVerbose(paste0("ID",journal.id,"_",title))
          }
          df1 <- data.frame(
            journal.id=real(journal.id),
            journals=real(title),
            abbr = real(subtitle),
            issn = real(issn),
            IF = real(IF),
            IF5 = real(IF5),
            selfcited = real(selfuse),
            h.index = real(h_index),
            CiteScore = real(score),
            disciplinary.quartile = real(q4),
            CiteScore.rank = real(rank),
            superior = real(superior),
            inferior = real(inferior),
            website = real(http),
            OA = real(oa),
            research = real(r),
            area  = real(area) ,
            regular = real(regular),
            BasicSubject = real(category.status[1]),
            SecondSubject = real(category.status[2]),
            TopJournal = real(category.status[3]),
            ReviewJournal = real(category.status[4]),
            speed = real(speed),
            recieve = real(recieve)
          )
        }

      } else {
        LuckyVerbose(a[["message"]])
        df1 <- data.frame(
          journal.id=real(journal.id),
          journals="Error",
          abbr = "Error",
          issn = "Error",
          IF = "Error",
          IF5 = "Error",
          selfcited = "Error",
          h.index = "Error",
          CiteScore = "Error",
          disciplinary.quartile = "Error",
          CiteScore.rank = "Error",
          superior = "Error",
          inferior = "Error",
          website = "Error",
          OA = "Error",
          research = "Error",
          area  = "Error" ,
          regular = "Error",
          BasicSubject = "Error",
          SecondSubject = "Error",
          TopJournal = "Error",
          ReviewJournal = "Error",
          speed = "Error",
          recieve = "Error"
        )

      }
    } else {
      df1 <- data.frame(
        journal.id=real(journal.id),
        journals="Error",
        abbr = "Error",
        issn = "Error",
        IF = "Error",
        IF5 = "Error",
        selfcited = "Error",
        h.index = "Error",
        CiteScore = "Error",
        disciplinary.quartile = "Error",
        CiteScore.rank = "Error",
        superior = "Error",
        inferior = "Error",
        website = "Error",
        OA = "Error",
        research = "Error",
        area  = "Error" ,
        regular = "Error",
        BasicSubject = "Error",
        SecondSubject = "Error",
        TopJournal = "Error",
        ReviewJournal = "Error",
        speed = "Error",
        recieve = "Error"
      )
      LuckyVerbose("Error:Not Available source!")
    }
    return(df1)
  }
  ano <- function(x)ano1(x,source,report)

  ##多个id同时提取
  df1 <- adply(as.matrix(journal.id),1,ano)
  df1 <- df1[,-1]
  colnames(df1)[match("journal.id",colnames(df1))] <- paste0(source,"id")

  ##输出结果
  if(report == T){LuckyVerbose("All done!")}
  return(df1)
}


####=========================MakeENlist===========================####
# MakeENlist create ENote journal list .txt to help annotation in ENote software.


#' @title create ENote journal list .txt to help annotation in ENote software
#' @description MakeENlist create ENote journal list .txt to help annotation in ENote software.
#' @keywords MakeENlist
#' @param data data frame or a matrix
#' @param journal.name.col the colname of journal full name
#' @param abbr.col the colname of journal abbreviation names
#' @param ISSN.col the colname of ISSN code
#' @param IF.col the colname of IF
#' @param names  part of saved file name
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## Use a journals data
#' MakeENlist(data = journals,
#'           journal.name.col="journals",
#'           abbr.col = "abbr",
#'           ISSN.col="issn",
#'           IF.col = "IF",
#'           IF.interval = "-",
#'           names = "love")
#' @export
MakeENlist <- function(data = journals,
                       journal.name.col="journals",
                       abbr.col = "abbr",
                       ISSN.col="issn",
                       IF.col = "IF",
                       IF.interval = "-",
                       names = "love"){
  ## 矩阵化
  data <- as.matrix(data)

  ## IF interval
  IF1 <- data[,IF.col]
  IF1 <- as.character(IF1)
  IF1 <- gsub("\\.",IF.interval,IF1)

  ## abbr首字母大写，其余小写
  #abbr1 = "abc"; convert.abbr(abbr1)
  #abbr1 = "abc def"; convert.abbr(abbr1)
  #abbr1 = "MACROMOL RAPID COMM"
  convert.abbr <- function(abbr1){
    a1 <- Fastextra(abbr1,"")
    p1 <- Fastgrep(c(" ","-"),a1)
    p1 <- p1[!is.na(p1)]
    if(length(p1)==0){
      #一个单词，无空格或-连接
      a1[1] <- toupper(a1[1])
      a1[2:length(a1)] <- tolower(a1[2:length(a1)])
      a2 <- paste(a1,collapse = "")
    } else {
      #多个单词有连接
      a1[c(1,p1+1)] <- toupper(a1[c(1,p1+1)])
      a1[setdiff(1:length(a1),c(1,p1+1))] <- tolower(a1[setdiff(1:length(a1),c(1,p1+1))])
      a2 <- paste(a1,collapse = "")
    }
    return(a2)
  }
  abbr <- apply(as.matrix(data[,abbr.col]),1,convert.abbr)

  ## multiple names
  fu.abbr <- abbr;
  fullnames <- as.character(data[,journal.name.col])
  fullnames1 <- paste(fu.abbr,as.character(data[,journal.name.col]),sep = ";  ")
  fn.if <-paste(fullnames1,IF1,sep = ";  IF:")
  fn.issn <- paste(fn.if,data[,ISSN.col],sep = ";  ISSN:")

  ## output
  df <- cbind(fullnames,fu.abbr,fn.if,fn.issn)
  write.table(df,paste0(names,"_Endnote.Journals.List.txt"),sep = "\t",row.names = F,col.names = F,quote = F)
  LuckyVerbose("Files have been saved!")

  ## End
}





