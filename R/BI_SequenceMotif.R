




#' @title co-motif exploration for sequence
#' @description co-motif exploration for sequence
#' @param sequence a vector of sequences
#' @param motif a charactor or numeric motif
#' @param base the type of base.Default is \code{c("A","C","U","G")}
#' @param verbose whether do reports
#' @importFrom plyr alply
#' @return a list contain counts and sequences information of motifs
#' @seealso \code{\link{combn2}};\code{\link[plyr]{alply}};
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## Create some random sequences
#' sequence <-  NULL
#' for(i in 1:100){
#'   set.seed(i);
#'   j <- sample(18:21,1)
#'   s.i <- sample(c("A","C","U","G"),j,replace = T)
#'   s.i <- paste(s.i,collapse = "")
#'   sequence <- c(sequence,s.i)
#' }
#'
#' ## explore two 4-base motifs
#' motif <- c("CCGU","GGAU")
#' l <- SequenceMotif(sequence,motif)
#' View(l$Counts)
#' View(l$Sequence)
#' l$Sequence$GGAU
#'
#' ## explore all 3-base or 4-base motifs
#' motif <- c(3,4)
#' l <- SequenceMotif(sequence,motif)
#' View(l$Counts)
#' View(l$Sequence)
#' l$Sequence$GGAU
#' @export
SequenceMotif <- function(sequence,
                          motif,
                          base = c("A","C","U","G"),
                          verbose = T){
  ## package
  #nd <- c("plyr");Plus.library(nd)

  ## get motif information
  sequence <- unique(sequence)
  if(!is.numeric(motif)){
    ## character motif
    LuckyVerbose("Character motif...")
    m1 <- sapply(motif, grepl,sequence)
  } else {
    ## numeric character
    if(verbose == T) LuckyVerbose("numeric motif...")
    motif2 <- combn2(base,motif)
    m1 <- sapply(motif2, grepl,sequence)
  }
  rownames(m1) <- sequence

  ## statistics
  m2 <- colSums(m1)
  m2 <- sort(m2,decreasing = T)
  m3 <- data.frame(motif = names(m2),counts = m2,portion = m2/length(sequence),stringsAsFactors = F)

  ## get
  if(verbose == T) LuckyVerbose("get sequences based on motifs...")
  m4 <- alply(.data = m1,
              .margins = 2,
              .fun = function(x) return(names(x)[x]))
  names(m4) <- colnames(m1)

  ## Output data
  l <- list(
    Counts = m3,
    Sequence = m4
  )
  if(verbose == T) LuckyVerbose("All done!")
  return(l)
}



#' @title get all arrage for given n and k
#' @description get all arrage for given n and k
#' @param vt a vector
#' @param n selected counts
#' @return a vector of arrage
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' base = c("A","C","U","G")
#' base2 <- combn2(base,c(3:5,9))
#' base3 <- combn2(base,c(3:5))
#' all(base2 %in% base3)
#' @export
combn2 <- function(vt,n){
  A <- NULL;a <- NULL
  for(i in n){
   if(i==1) for (b1 in vt) a <- c(a,paste0(b1,collapse =""))
   if(i==2) for (b1 in vt) for (b2 in vt) a <- c(a,paste0(b1,b2,collapse =""))
   if(i==3) for (b1 in vt) for (b2 in vt) for (b3 in vt) a <- c(a,paste0(b1,b2,b3,collapse =""))
   if(i==4) for (b1 in vt) for (b2 in vt) for (b3 in vt) for (b4 in vt) a <- c(a,paste0(b1,b2,b3,b4,collapse =""))
   if(i==5) for (b1 in vt) for (b2 in vt) for (b3 in vt) for (b4 in vt) for (b5 in vt) a <- c(a,paste0(b1,b2,b3,b4,b5,collapse =""))
   if(i==6) for (b1 in vt) for (b2 in vt) for (b3 in vt) for (b4 in vt) for (b5 in vt) for (b6 in vt) a <- c(a,paste0(b1,b2,b3,b4,b5,b6,collapse =""))
   if(i==7) for (b1 in vt) for (b2 in vt) for (b3 in vt) for (b4 in vt) for (b5 in vt) for (b6 in vt) for (b7 in vt) a <- c(a,paste0(b1,b2,b3,b4,b5,b6,b7,collapse =""))
   if(i==8) for (b1 in vt) for (b2 in vt) for (b3 in vt) for (b4 in vt) for (b5 in vt) for (b6 in vt) for (b7 in vt) for (b8 in vt) a <- c(a,paste0(b1,b2,b3,b4,b5,b6,b7,b8,collapse =""))
   if(!i %in% 1:8) LuckyVerbose(paste0(i," is not in 1 to 8,which is not supported."))

   ## 汇总
   A <- c(A,a)
   a <- NULL
  }
  return(A)
}


















