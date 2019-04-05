


#' @title Give Mutiple primers,test the amplication sequence and whether cross different exons
#' @description Give Mutiple primers,test the amplication sequence and whether cross different exons
#' @param forward a name vector of forward primers
#' @param reverse a name vector of reverse primers
#' @param path.transcript.fa the path of transcript(cDNA) fa file.if primers are multiple,it is recommanded to set \code{path.transcript.fa=NULL}
#' @param path.exon.fa the path of exons fa file.if primers are multiple,it is recommanded to set \code{path.exon.fa=NULL}
#' @importFrom Biostrings readDNAStringSet reverseComplement subseq vmatchPattern DNAString oligonucleotideFrequency
#' @importFrom httr GET content_type stop_for_status content
#' @return A luckyList
#' @seealso \code{\link[Biostrings]{XStringSet-io}};\code{\link[XVector]{XVector-class}};\code{\link[Biostrings]{matchPattern}};\code{\link{testPrimer1}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' forward = c(DSCR8="GGGTGGCAAAAAGAGAAGATGG",
#'             DSCR8 = "CGTCTTGTCCACTCGTTTCTG",
#'             DSCR8 = "AAATAGGCTTCCCACTTGGC",
#'             DSCR8 = "GGGTGGCAAAAAGAGAAGATGG")
#' reverse = c(DSCR8="GGTCCAGGCTCCTTCATTTTTG",
#'             DSCR8 = "GGCATGTTCTAAATGGCAGCTC",
#'             DSCR8 = "CCAGGCTCCTTCATTTTTGCTC",
#'             DSCR8 = "GGTCCAGGCTCCTTCATTTTTG")
#' res <- testPrimer(forward,reverse)
#'
#' library(lucky)
#' forward = c(FASN = "ACAGCGGGGAATGGGTACT",
#'             SCD = "TCATAATTCCCGACGTGGCT",
#'             PPARG = "GGGATCAGCTCCGTGGATCT",
#'             ACACA = "TCACACCTGAAGACCTTAAAGCC",
#'             ACTB = "CATGTACGTTGCTATCCAGGC",
#'             TNF = "TGGCGTGGAGCTGAGAGATA",
#'             IL6 = "GAGTAGTGAGGAACAAGCCAGA",
#'             PCK1 = "AAGGTGGGGACGTTCAGAATC",
#'             G6PC = "CACTTCCGTGCCCCTGATAA")
#' reverse = c(FASN = "GACTGGTACAACGAGCGGAT",
#'             SCD = "CCCAGAAATACCAGGGCACA",
#'             PPARG = "TGCACTTTGGTACTCTTGAAGTT",
#'             ACACA = "AGCCCACACTGCTTGTACTG",
#'             ACTB = "CTCCTTAATGTCACGCACGAT",
#'             TNF = "TGATGGCAGAGAGGAGGTTG",
#'             IL6 = "AAGCTGCGCAGAATGAGATGA",
#'             PCK1 = "ATGGTTCCCTGCCTCATTCC",
#'             G6PC = "TAGTATACACCTGCTGTGCCC")
#' res <- testPrimer(forward,reverse)
#' @export
testPrimer <- function(forward,
                       reverse,
                       path.transcript.fa=NULL,
                       path.exon.fa=NULL,
                       verbose = T){
  ## load package
  #nd <- c("XML","rvest","httr","jsonlite","xml2","Biostrings")
  #Plus.library(nd)

  ## gene.id
  gene <- unique(names(forward))
  gene2 <- NULL
  for(i in 1:length(gene)){
    g.i <- gene[i]
    if(length(grep("ENSG",g.i)) == 0){
      g.i <- convert(g.i,fromtype = "SYMBOL",totype = "ENSEMBL")
    } else {
      g.i <- g.i
    }
    gene2 <- c(gene2,g.i)
  }
  gene <- gene2

  ## Create a new file to place .fa
  old <- options()
  options(warn = -1)
  dir.create("testPrimer",recursive = T)
  on.exit(options(old), add = TRUE)
  path1 <- list.files(path = "./testPrimer",pattern = ".fa",full.names = T)
  path2 <- NULL
  for(i in 1:length(gene)){ # i=1
    path2.i <- list.files(path = "./testPrimer",pattern = gene[i],full.names = T)
    path2 <- c(path2,path2.i)
  }
  path <- intersect(path1,path2)
  if(length(path)==0){

    ## Get .fa and path
    if(verbose) LuckyVerbose("Download FASTA from http://rest.ensembl.org...")

    if(is.null(path.transcript.fa)|is.null(path.exon.fa)){
      ## get transcrips and exon names of given gene.
      transcrips <- NULL;exons <- NULL
      for(i in 1:length(gene)){ # i=1
        g.i <- gene[i]
        tran.i <- common.annot.transcripts_GRCh38$transcript_id[common.annot.transcripts_GRCh38$gene_id %in% g.i]
        tran.i <- unique(tran.i[!is.na(tran.i)])
        exon.i <-  common.annot.transcripts_GRCh38$exon_id[common.annot.transcripts_GRCh38$gene_id%in% g.i]
        exon.i <- unique(exon.i[!is.na(exon.i)])
        transcrips <- c(transcrips,list(tran.i));names(transcrips)[i] <- g.i
        exons <- c(exons,list(exon.i));names(exons)[i] <- g.i
      }

      ## Download fasta files via API from ensembl
      server <- "http://rest.ensembl.org"
      for(i in 1:length(gene)){ # i=1
        g.i <- gene[i]
        tran.i <- transcrips[[g.i]]
        tran.fa <- NULL;
        for(j in 1:length(tran.i)){ # j=1
          tran.j <- tran.i[j]
          # ext <- "/sequence/id/ENST00000288602?type=cdna"
          ext <- paste0("/sequence/id/",tran.j,"?type=cdna")
          r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("text/x-fasta"))
          httr::stop_for_status(r)
          fa.j <- httr::content(r,as = "text")
          tran.fa <- c(tran.fa,fa.j)
        }
        write.table(tran.fa,file = paste0("./testPrimer/",g.i,"_",convert(g.i),"_transcrips.fa"),quote = F,col.names = F,row.names = F)
        ## exons
        exon.i <- exons[[i]]
        exon.fa <- NULL
        for(j in 1:length(exon.i)){ # j=1
          exon.j <- exon.i[j]
          # ext <- "/sequence/id/ENSE00001154485?type=genomic"
          ext <-  paste0("/sequence/id/",exon.j,"?type=genomic")
          r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("text/x-fasta"))
          httr::stop_for_status(r)
          fa.j <- httr::content(r,as = "text")
          exon.fa <- c(exon.fa,fa.j)
        }
        write.table(exon.fa,file = paste0("./testPrimer/",g.i,"_",convert(g.i),"_exons.fa"),quote = F,col.names = F,row.names = F)

      }

      ## create path for testPrimer1
      path1 <- list.files(path = "./testPrimer",pattern = ".fa",full.names = T)
      path2 <- NULL
      for(i in 1:length(gene)){ # i=1
        path2.i <- list.files(path = "./testPrimer",pattern = gene[i],full.names = T)
        path2 <- c(path2,path2.i)
      }
      path <- intersect(path1,path2)
      tran.path <- path[grep("transcrips",path)]
      exon.path <- path[grep("exons",path)]

    } else {
      tran.path <- path.transcript.fa
      exon.path <- path.exon.fa
    }
    if(verbose) LuckyVerbose('FASTA had been saved in "./testPrimer"!')
  } else {
    if(verbose) LuckyVerbose("Located FASTA exited and would be loaded!")
    tran.path <- path[grep("transcrips",path)]
    exon.path <- path[grep("exons",path)]
  }

  ## test primer
  l <- NULL
  for(j in 1:length(forward)){ # j=1
    forward.i <- as.character(forward[j])
    reverse.i <- as.character(reverse[j])

    # gene name
    g.i <- names(forward[j])
    if(length(grep("ENSG",g.i)) == 0){
      g.i <- convert(g.i,fromtype = "SYMBOL",totype = "ENSEMBL")
    } else {
      g.i <- g.i
    }

    # path
    path.p <- tran.path[grep(g.i,tran.path)]
    exon.p <- exon.path[grep(g.i,exon.path)]

    # testPrimer1
    l.i <- testPrimer1(forward=forward.i,
                       reverse=reverse.i,
                       path.transcript.fa = path.p,
                       path.exon.fa = exon.p,
                       verbose = F)
    n.i <- paste(g.i,convert(g.i),paste0("F:",forward.i),paste0("R:",reverse.i),sep = "_")
    am.i <- l.i[["Data"]][["UniqueAmplication"]]

    # report
    if(l.i$Data$CrossExons == "Cross different exons"){
      report.exon <- "Cross different exons"
    } else {
      report.exon <- NULL
    }
    if(verbose) LuckyVerbose(paste0(paste0(j,".",n.i),": There are ",length(am.i)," sequence amplicated by this primer pair;",report.exon),levels = 1)
    l <- c(l,list(l.i));names(l)[j] <- n.i

  }
  return(l)
}

#' @title Give a primer,test the amplication sequence and whether cross different exons
#' @description Give a primer,test the amplication sequence and whether cross different exons
#' @param forward forward primer
#' @param reverse reverse primer
#' @param path.transcript.fa the path of transcript(cDNA) fa file
#' @param path.exon.fa the path of exons fa file
#' @importFrom Biostrings readDNAStringSet reverseComplement subseq vmatchPattern DNAString oligonucleotideFrequency
#' @importFrom httr GET
#' @return A luckyList
#' @seealso \code{\link[Biostrings]{XStringSet-io}};\code{\link[XVector]{XVector-class}};\code{\link[Biostrings]{matchPattern}};
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' forward = "AAGGAGCCTGGACCCAACTT"
#' reverse = "TGAGCCAGCCTTCATTTTCCA"
#' path.transcipt.fa = "E:/iProjects/XYJ.lnc/constant/Homo_sapiens_DSCR8_cDNA.fa"
#' path.exon.fa = "E:/iProjects/XYJ.lnc/constant/Homo_sapiens_DSCR8_Exons.fa"
#'
#' res <- testPrimer1(forward,
#'                    reverse,
#'                    path.transcript.fa,
#'                    path.exon.fa)
#' ## see result
#' names(res$Data$Amplication)
#' res$Data$CrossExons
#'
#' @export
testPrimer1 <- function(forward,
                        reverse,
                        path.transcript.fa,
                        path.exon.fa,
                        verbose = T){
  ## load package
  #nd <- c("Biostrings")
  #Plus.library(nd)

  ## Add cDNA sequences and exon sequences
  cDNA <- readDNAStringSet(path.transcript.fa, format="fasta")
  exon <- readDNAStringSet(path.exon.fa, format="fasta")

  ## match amplication sequence
  am_f <- as.list(vmatchPattern(forward,cDNA))
  am_r <- as.list(vmatchPattern(reverseComplement(DNAString(reverse)),cDNA)) # reverse complement sequence

  ## find out start and end;Output the amplication sequence
  get1 <- function(am_f.i,am_r.i,cDNA.i){
    # am_f.i <- am_f[[2]]
    # am_r.i <- am_r[[2]]
    # cDNA.i <- cDNA[[2]]
    start <- am_f.i@start
    end <- am_r.i@start + am_r.i@width - 1
    if(length(end) == 0 | length(start) == 0){
      ## no amplication sequence in this transcript
      amplication.i <- NA
    } else {
      amplication.i <- subseq(cDNA.i,start = start,end = end)
    }
    return(amplication.i)
  }
  amplication <- NULL
  for(i in 1:length(cDNA)){ # i=1
    am_f.i <- am_f[[i]]
    am_r.i <- am_r[[i]]
    cDNA.i <- cDNA[[i]]
    a.i <- get1(am_f.i,am_r.i,cDNA.i)
    amplication <- c(amplication,list(a.i))
    names(amplication)[i] <- names(cDNA)[i]
  }
  test <- unlist(amplication)
  amplication <- amplication[!is.na(test)]

  if(length(amplication) == 0){
    # no amplication
    if(verbose) LuckyVerbose(paste0("There are no sequence amplicated by this primer pair."),levels = 1)
  } else {
    ## get unique amplication
    am_2 <- NULL
    for(i in 1:length(amplication)){ # i=1
      a.i <- amplication[[i]]
      a.i <- as.character(unlist(a.i))
      am_2 <- c(am_2,a.i);
    }
    am_2 <- unique(am_2)
    if(verbose) LuckyVerbose(paste0("There are ",length(am_2)," sequence amplicated by this primer pair."),levels = 1)

    ## Calculate Tm value
    if(F){
      am_3 <- am_2
      Tm <- NULL
      for(i in 1:length(am_2)){ # i=1
        a.i <- am_3[i]
        freq.i <- oligonucleotideFrequency(DNAString(a.i),width = 1)
        CG <- freq.i[names(freq.i) %in% c("G","C")]
        CG.percentage <- sum(CG)/sum(freq.i)
        # 对于更长的寡聚核苷酸，Tm计算公式为：Tm = 81.5 + 16.6 x Log10[Na+] + 0.41 (%GC) – 600/size公式中，Size = 引物长度
        Tm.i <- CG.percentage/2.44 + 69.3
      }
    }

    ## whether cross different exons
    get2 <- function(am.i,exon){
      # am.i <- am_2[i]
      test.i <- as.list(vmatchPattern(pattern = am.i,subject = exon))
      cross <- NULL
      for(i in 1:length(exon)){ # i=1
        start <- test.i[[i]]@start
        if(!length(start) == 0){
          c.i <- test.i[[i]]
          cross <- c(cross,list(c.i))
          names(cross)[length(cross)] <- names(exon)[i]
        }
      }
      return(cross)
    }
    cross <- NULL
    for(i in 1:length(am_2)){ # i =1
      a.i <- am_2[i]
      cross.i <- get2(a.i,exon)
      if(is.null(cross.i)){
        ## cross different exons
        cross.i <- "Cross different exons"
      } else {
        cross.i <- cross.i
      }
      cross <- c(cross,list(cross.i));
      names(cross)[i] <- a.i
    }

    ## Output Data
    l <- list(
      Repeat = list(
        forward = forward,
        reverse = reverse,
        path.transcript.fa = path.transcript.fa,
        path.exon.fa = path.exon.fa
      ),
      Data = list(
        transcript = cDNA@ranges@NAMES,
        Amplication = amplication,
        UniqueAmplication =  am_2,
        CrossExons = cross
      )
    )
    return(l)
  }

}























