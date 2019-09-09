#' Generate a list of genomic features GRangesList objects
#' 
#' This function extracts and filters gene, upstream and downstream features from a \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object.
#' @param first A \code{\link[GenomicRanges:GenomicRanges-class]{GenomicRanges}} or \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} object.
#' @param second An optional \code{\link[GenomicRanges:GenomicRanges-class]{GenomicRanges}} or \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param method A character string indicating which reference should be used to assign the \code{\link{seqlevels}}. One of 'first' (default), 'second', 'union', "txdb": can be abbreviated.
#' @return A list of \code{\link[GenomicRanges:GenomicRanges-class]{GenomicRanges}} or \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} object(s).
#'
#' @import IRanges
#' @import BiocGenerics
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import AnnotationDbi
#' @export
#' @examples
#' library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
#' 
#' TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
#'
#' first = GRanges(c("chr4", "chr2L"), IRanges(c(10, 20), c(30, 30)), c("-", "+"))
#' second = GRanges(c("chr4", "chr3R"), IRanges(c(10, 20), c(30, 30)), c("-", "+"))
#'
#' # Assign seqlevels and seqinfo from a TxDb object
#' fixList = fixSeqInfo(first, TxDb=TxDb)
#' sapply(fixList, head)
#' sapply(fixList, seqlevels)
#' sapply(fixList, seqinfo)
#'
#' # Assign the seqlevels of the first to the second object
#' fixList = fixSeqInfo(first, second, method="first")
#' sapply(fixList, head)
#' sapply(fixList, seqlevels)
#' sapply(fixList, seqinfo)
#'
#' # Assign the union of the seqlevels to both objects
#' fixList = fixSeqInfo(first, second, method="union")
#' sapply(fixList, head)
#' sapply(fixList, seqlevels)
#' sapply(fixList, seqinfo)


fixSeqInfo <- function(first=NULL, second=NULL, TxDb=NULL, method=c("first", "second", "union", "txdb")){

   # Check for required args
   stopifnot( !is.null(first) & (!is.null(second) | !is.null(TxDb)) )
   method = match.arg(method)

   seqlevelsStyle(first) = "UCSC"
   first = sortSeqlevels(keepStandardChromosomes(first, pruning.mode="tidy"))
   if( !is.null(second) ){
      seqlevelsStyle(second) = "UCSC"
      second = sortSeqlevels(keepStandardChromosomes(second, pruning.mode="tidy"))
   }else if( method%in%c("second", "union") ){
      stop("The method selected requires a second object!")
   }
   if( !is.null(TxDb) ){
      seqlevelsStyle(TxDb) = "UCSC"
      TxDb = sortSeqlevels(keepStandardChromosomes(TxDb, pruning.mode="tidy"))
      if( !all(seqlevels(first)%in%seqlevels(TxDb)) ){
         stop("Some seqlevels of the first object are not present in the TxDb object!")
      }
      if( !is.null(second) ){
         if( !all(seqlevels(second)%in%seqlevels(TxDb)) ){
            stop("Some seqlevels of the second object are not present in the TxDb object!")
         }
      }
   }else if( method == "txdb" ){
         stop("TxDb object is missing!")
   }

   if( !is.null(second) ){
      if( method == "first" ){
         extraLevels = setdiff(seqlevels(second), seqlevels(first))
         second = dropSeqlevels(second, extraLevels, pruning.mode="tidy")
         seqlevels(second) = seqlevels(first)
      }else if( method == "second" ){
         extraLevels = setdiff(seqlevels(first), seqlevels(second))
         first = dropSeqlevels(first, extraLevels, pruning.mode="tidy")
         seqlevels(first) = seqlevels(second)
      }else if( method == "union" ){
         sharedLevels = union(seqlevels(first), seqlevels(second))
         seqlevels(first) = sharedLevels
         seqlevels(second) = sharedLevels
      }
      if( !is.null(TxDb) ){
         if( method == "txdb" ){
            seqlevels(first) = seqlevels(TxDb)
            first = keepSeqlevels(first, seqlevels(TxDb), pruning.mode="tidy")
            seqlevels(second) = seqlevels(TxDb)
            second = keepSeqlevels(second, seqlevels(TxDb), pruning.mode="tidy")
         }
         seqinfo(first) = seqinfo(TxDb)[seqlevels(first)]
         seqinfo(second) = seqinfo(TxDb)[seqlevels(second)]
      }
      return(list(first, second))
   }else{
      if( method == "txdb" ){
         seqlevels(first) = seqlevels(TxDb)
         first = keepSeqlevels(first, seqlevels(TxDb), pruning.mode="tidy")
      }
      seqinfo(first) = seqinfo(TxDb)[seqlevels(first)]
      return(list(first))
   }
}
