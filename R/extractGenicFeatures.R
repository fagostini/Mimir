#' Generate a list of genic features GRangesList objects
#' 
#' This function extracts and filters the genic features from a \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param tx2gene A \code{"\linkS4class{data.frame}"} object. The first column must be of transcript identifiers, while the second must be of gene identifiers. Additional columns will be discarded.
#' @param selectGn A vector of optional gene identifiers to keep.
#' @param selectTx A vector of optional transcript identifiers to keep.
#' @param excludeIntrons When set to ‘TRUE’, the extraction of intronic regions is skipped.
#' @param bins A 4 integers ordered vector. The vector order determines the 5'UTR, CDS, 3'UTR and Introns minumum region widths.
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
#' genicRegions = extractGenicFeatures(TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene)

extractGenicFeatures <- function(TxDb=NULL, tx2gene=NULL, selectGn=NULL, selectTx=NULL, excludeIntrons=TRUE, bins=c(20, 100, 70, 100), verbose=TRUE){

    # Check for required args
    stopifnot( !is.null(TxDb) & length(bins)==4 )
    seqlevelsStyle(TxDb) = "UCSC"

    # Generate transcript -> gene reference table
    if( is.null(tx2gene) ){
        if( verbose ) message("Generating Reference Table...")
        tx2gene = unlist(transcriptsBy(TxDb, by="gene"), use.names=TRUE)
        tx2gene = data.frame(tx_id = tx2gene$tx_name, gene_id = names(tx2gene))
    }else{
        tx2gene = data.frame(tx2gene)
        tx2gene = tx2gene[, 1:2]
        names(tx2gene) = c("tx_id", "gene_id")
    }

    # Select by gene
    if( !is.null(selectGn) ){
        tx2gene = tx2gene[gene_id%in%selectGn,]
    }
    # Select by transcript
    if( !is.null(selectTx) ){
        tx2gene = tx2gene[tx_id%in%selectTx,]
    }

    if( verbose )
        message(paste("Using",
            length(unique(tx2gene$gene_id)), "Genes and",
            length(unique(tx2gene$tx_id)), "Transcripts"))

    # Extract 5'UTRs
    fiveUTRs = fiveUTRsByTranscript(TxDb, use.names=TRUE)
    fiveUTRs = fiveUTRs[names(fiveUTRs)%in%tx2gene$tx_id]
    fiveUTRs = fiveUTRs[sum(width(fiveUTRs))>=bins[1]]
    seqinfo(fiveUTRs) = seqinfo(TxDb)
    fiveUTRs = keepStandardChromosomes(fiveUTRs, pruning.mode="coarse")
    if( verbose ) message(paste("Extracted", length(fiveUTRs), "5'UTRs"))

    # Extract coding sequences
    cds = cdsBy(TxDb, by="tx", use.names=TRUE)
    cds = cds[names(cds)%in%tx2gene$tx_id]
    cds = cds[sum(width(cds))>=bins[2]]
    seqinfo(cds) = seqinfo(TxDb)
    cds = keepStandardChromosomes(cds, pruning.mode="coarse")
    if( verbose ) message(paste("Extracted", length(cds), "CDSs"))

    # Extract 3'UTRs
    threeUTRs = threeUTRsByTranscript(TxDb, use.names=TRUE)
    threeUTRs = threeUTRs[names(threeUTRs)%in%tx2gene$tx_id]
    threeUTRs = threeUTRs[sum(width(threeUTRs))>=bins[3]]
    seqinfo(threeUTRs) = seqinfo(TxDb)
    threeUTRs = keepStandardChromosomes(threeUTRs, pruning.mode="coarse")
    if( verbose ) message(paste("Extracted", length(threeUTRs), "5'UTRs"))

    genicRegions = list(UTR5 = fiveUTRs, CDS = cds, UTR3 = threeUTRs)

    if( !excludeIntrons ){
        # Extract intronic sequences
        if( verbose ) message("Extracting Introns...")
        introns = intronsByTranscript(TxDb, use.names=TRUE) 
        introns = introns[names(introns)%in%tx2gene$tx_id]
        introns = introns[sum(width(introns))>=bins[4]]
        seqinfo(introns) = seqinfo(TxDb)
        introns = keepStandardChromosomes(introns, pruning.mode="coarse")
        if( verbose ) message(paste("Extracted", length(introns), "introns"))

        genicRegions = append(genicRegions, list(Intron = introns))
    }

    return(genicRegions)
}
