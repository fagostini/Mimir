#' Generate a list of genomic features GRangesList objects
#' 
#' This function extracts and filters gene, upstream and downstream features from a \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param selectGn A vector of optional gene identifiers to keep.
#' @param excludeIntrons When set to 'TRUE', the extraction of intronic regions is skipped.
#' @param exon_width A positive integer. It determines the minumum width for the sum of all exons in a gene.
#' @param intron_width A positive integer. It determines the minumum width for the sum of all introns in a gene.
#' @param upstream_width A positive integer. It determines the upstream region width.
#' @param downstream_width A positive integer. It determines the downstream gene body width.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
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

extractGenomicFeatures <- function(TxDb=NULL, selectGn=NULL, excludeIntrons=TRUE,
    exon_width=1000, intron_width=1000, upstream_width=1000, downstream_width=1000, verbose=TRUE){

    # Check for required args
    stopifnot( !is.null(TxDb) )
    seqlevelsStyle(TxDb) = "UCSC"


    # Extract either whole genes or gene exons
    if( !excludeIntrons ){
        exons = genes(TxDb, columns="gene_id")
        names(exons) = NULL
        exons = reduce(split(exons, exons$gene_id))
    }else{
       exons = reduce(exonsBy(TxDb, by="gene"))
    }

   exons = exons[sum(width(exons))>=exon_width]
   seqinfo(exons) = seqinfo(TxDb)
   exons = keepStandardChromosomes(exons, pruning.mode="coarse")

    # Select by gene
    if( !is.null(selectGn) ){
         exons = exons[names(exons)%in%selectGn]
    }
    if( verbose ) message(paste("Extracted", length(exons), "genes"))

    # Extract upstream regions
    upstreams = genes(TxDb, columns="gene_id", filter=list(gene_id=names(exons)))
    names(upstreams) = NULL
    upstreams = reduce(split(upstreams, upstreams$gene_id))
    upstreams = flank(upstreams, width=upstream_width, start=TRUE, both=FALSE, use.names=TRUE)
    seqinfo(upstreams) = seqinfo(TxDb)
    upstreams = subsetByOverlaps(upstreams, as(seqinfo(TxDb), "GRanges"), type="within")
    if( verbose ) message(paste("Extracted", length(upstreams), "upstream regions"))

    # Extract downstream regions
    downstreams = genes(TxDb, columns="gene_id", filter=list(gene_id=names(exons)))
    names(downstreams) = NULL
    downstreams = reduce(split(downstreams, downstreams$gene_id))
    downstreams = flank(downstreams, width=downstream_width, start=FALSE, both=FALSE, use.names=TRUE)
    seqinfo(downstreams) = seqinfo(TxDb)
    downstreams = subsetByOverlaps(downstreams, as(seqinfo(TxDb), "GRanges"), type="within")
    if( verbose ) message(paste("Extracted", length(downstreams), "downstream regions"))

    genomicRegions = list(Upstream = upstreams,
                          Exon = exons,
                          Downstream = downstreams)

    if( !excludeIntrons ){
        introns = genes(TxDb)
        introns = split(introns, names(introns))
        introns = reduce(introns[names(exons)])
        introns = setdiff(introns, exons)
        introns = introns[sum(width(introns))>=intron_width]
        genomicRegions = append(genomicRegions, list(Intron = introns))
    }

    return(genomicRegions)
}
