#' Generate a list of genomic features GRangesList objects
#' 
#' This function extracts and filters gene, upstream and downstream features from a \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param selectGn A vector of optional gene identifiers to keep.
#' @param excludeIntrons When set to ‘TRUE’, the extraction of intronic regions is skipped.
#' @param body_width A positive integer. It determines the minumum gene body width.
#' @param upstream_width A positive integer. It determines the upstream region width.
#' @param downstream_width A positive integer. It determines the downstream gene body width.
#' @param verbose When set to ‘TRUE’, the function prints diagnostic messages.
#' @return A named list of \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} objects.
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
#' genomicRegions = extractGenomicFeatures(TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene)

extractGenomicFeatures <- function(TxDb=NULL, selectGn=NULL, excludeIntrons=TRUE, body_width=1000, upstream_width=100, downstream_width=1000, verbose=TRUE){

    # Check for required args
    stopifnot( !is.null(TxDb) )
    seqlevelsStyle(TxDb) = "UCSC"


    # Extract either whole genes or gene exons
    if( !excludeIntrons ){
        bodies = genes(TxDb, columns="gene_id")
        names(bodies) = NULL
        bodies = reduce(split(bodies, bodies$gene_id))
    }else{
       bodies = reduce(exonsBy(TxDb, by="gene"))
    }

   bodies = bodies[sum(width(bodies))>=body_width]
   seqinfo(bodies) = seqinfo(TxDb)
   bodies = keepStandardChromosomes(bodies, pruning.mode="coarse")

    # Select by gene
    if( !is.null(selectGn) ){
         bodies = bodies[names(bodies)%in%selectGn]
    }
    if( verbose ) message(paste("Extracted", length(bodies), "genes"))

    # Extract upstream regions
    upstreams = genes(TxDb, columns="gene_id", filter=list(gene_id=names(bodies)))
    names(upstreams) = NULL
    upstreams = reduce(split(upstreams, upstreams$gene_id))
    upstreams = flank(upstreams, width=upstream_width, start=TRUE, both=FALSE, use.names=TRUE)
    seqinfo(upstreams) = seqinfo(TxDb)
    upstreams = keepStandardChromosomes(upstreams, pruning.mode="coarse")
    if( verbose ) message(paste("Extracted", length(upstreams), "upstream regions"))

    # Extract downstream regions
    downstreams = genes(TxDb, columns="gene_id", filter=list(gene_id=names(bodies)))
    names(downstreams) = NULL
    downstreams = reduce(split(downstreams, downstreams$gene_id))
    downstreams = flank(downstreams, width=downstream_width, start=FALSE, both=FALSE, use.names=TRUE)
    seqinfo(downstreams) = seqinfo(TxDb)
    downstreams = keepStandardChromosomes(downstreams, pruning.mode="coarse")
    if( verbose ) message(paste("Extracted", length(downstreams), "downstream regions"))

    genomicRegions = list(Upstream = upstreams, Gene_body = bodies, Downstream = downstreams)

    return(genomicRegions)
}
