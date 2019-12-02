#' Generate a GRanges object with exon/intron classification
#' 
#' This function extracts and labels the exonic (and intronic) features according to their rank within the genes/transcripts.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param by A character string indicating whether features should be extracted at gene ('gene') or transcript ('tx') level.
#' @param selectGn A vector of optional gene/transcript identifiers to keep.
#' @param excludeIntrons When set to 'TRUE', the extraction of intronic regions is skipped.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
#' @return A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} object.
#'
#' @import IRanges
#' @import BiocGenerics
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import AnnotationDbi
#' @import data.table
#' @export
#' @examples
#' library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
#' 
#' genomicRegions = extractGeneFeatureGroups(TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene)

extractGeneFeatureGroups <- function(TxDb=NULL, by=c("gene", "tx"), selectGn=NULL, excludeIntrons=TRUE, verbose=TRUE){

    # Check for required args
    stopifnot( !is.null(TxDb) )
    seqlevelsStyle(TxDb) = "UCSC"

    by = match.arg(by)

    # Extract exons only or exons and introns
    if( excludeIntrons ){
        exons = reduce(exonsBy(TxDb, by=by))
        exons = exons[sapply(exons, length)>2]
    }else{
       exons = reduce(exonsBy(TxDb, by=by))
       exons = exons[sapply(exons, length)>2]
       if( by == "gene"){
            introns = genes(TxDb, columns="gene_id")
            names(introns) = NULL
            introns = reduce(split(introns, introns$gene_id))
       }else{
            introns = transcriptsBy(TxDb, by="gene")
       }
       introns = reduce(setdiff(introns[names(exons)], exons))
    }
    
    # Select by gene/transcript id
    if( !is.null(selectGn) ){
         exons = exons[names(exons)%in%selectGn]
         introns = introns[names(introns)%in%selectGn]
    }
    if( verbose ) message(paste("Extracted", length(exons), "transcripts"))
    if( verbose ) message(paste("    exons:", length(unlist(exons))))
    if( verbose & !excludeIntrons ) message(paste("  introns:", length(unlist(introns))))
    if( verbose ) warning("Extracted features have not been filtered by size!")

    if( excludeIntrons ){
        exons = list(Exon=exons)
    }else{
        exons = list(Exon=exons, Intron=introns)
    }

    exons = rbindlist(lapply(exons,
        function(x)
            rbindlist(lapply(x,
                function(f)
                    data.table(as.data.frame(f))),
                idcol=paste(by, "id", sep="_"))),
            idcol="feature")

    if( !excludeIntrons )
        exons[feature=="Intron", label := "Introns"]
    setkeyv(exons, c(paste(by, "id", sep="_"), "seqnames", "start", "end"))
    setkeyv(exons, paste(by, "id", sep="_"))
    
    assignLabel <- function(v, n, s){
        res = sapply(v,
            function(i){
                if( i == 1 ){
                    if( s == "+" ){
                        return("First_exon")
                    }else{
                        return("Last_exon")
                    }
                }else if( i == n ){
                    if( s == "+" ){
                        return("Last_exon")
                    }else{
                        return("First_exon")
                    }
                }else{
                    return("Middle_exons")
                }
            }
        )
        return(res)
    }

    exons[feature=="Exon", label := assignLabel(sequence(.N), .N, unique(strand)), by=key(exons)]

    exons = with(exons, GRanges(seqnames, IRanges(start, end), strand, eval(parse(text = paste(by, "id", sep="_"))), label))

    colnames(mcols(exons)) = c(paste(by, "id", sep="_"), "label")

    return(exons)
}
