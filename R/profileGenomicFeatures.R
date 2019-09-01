#' Calculate metadata profiles across genomic features
#' 
#' This function calculates metadata profiles across genomic features using a procedure similar to that used in \href{https://www.sciencedirect.com/science/article/pii/S1097276519303533?via%3Dihub#fig1}{Viphakone et al., 2019}.
#' @param genomicRegions A named list of \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} objects.
#' @param sampleObject A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} or \code{\link[GenomicAlignments:GAlignments-class]{GenomicAlignments}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. Required if \code{genomicRegions} is not provided. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param bins An ordered integer vector (must be greater equal that the length of \code{genomicRegions}), or 3 if the latter is not provided. The vector order determines the number of bins for each region. If more bins than regions are provided, the additional will be ignored.
#' @param weightCol A single character string. This must be the name of an integer column in the \code{sampleObject} object.
#' @param ignoreStrand When set to 'TRUE', the strand information in \code{sampleObject} is ignored. This does not affect the features in \code{genomicRegions}.
#' @param dropEmpty When set to 'TRUE', the transcripts with no signal in any of their sub-regions will be discarded. When set to 'FALSE', all values of these regions will be set to 0.
#' @param collapse When set to 'TRUE', the profiles are collapsed into a single profile.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
#' @return A \code{\link[data.table:data.table-class]{data.table}} of the normalised binned coverage across the genomic features. Column names are determined by \code{collapse}.
#'
#' @import IRanges
#' @import BiocGenerics
#' @importFrom S4Vectors mcols
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom data.table ":=" ".N" ".GRP" ".I" data.table as.data.table setkey setkeyv setnames rbindlist
#' @importFrom stats setNames
#' @export
#' @examples
#' library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
#' TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
#' 
#' query = extractGenomicFeatures(TxDb)
#' 
#' library("pasillaBamSubset")
#' library("GenomicAlignments")
#' 
#' fl1 <- untreated1_chr4()
#' subject = readGAlignments(fl1)
#'
#' profile = profileGenomicFeatures(genomicRegions=query, sampleObject=subject, TxDb=TxDb)
#'
#' library("ggplot2")
#'
#' ggplot(profile, aes(x=bin, y=Mean, col=region_id)) + 
#'    geom_ribbon(aes(ymin=Mean-Sd, ymax=Mean+Sd, fill=region_id), alpha=0.5) +
#'    geom_line() +
#'    geom_vline(xintercept=c(10.5, 110.5), linetype="dashed") +
#'    scale_x_continuous("Relative position") +
#'    scale_y_continuous("Average normalised signal")

profileGenomicFeatures <- function(genomicRegions=NULL, sampleObject=NULL, bins=c(10, 100, 100), TxDb=NULL, weightCol=NULL, ignoreStrand=FALSE, dropEmpty=TRUE, collapse=TRUE, verbose=TRUE){

        stopifnot( !is.null(genomicRegions) | !is.null(TxDb) )
        stopifnot( !is.null(genomicRegions) & length(genomicRegions)<=length(bins) ) # Early call
        stopifnot( !is.null(weightCol)  & weightCol%in%colnames(mcols(sampleObject)) )
        stopifnot( !any(is.na(seqlengths(sampleObject))) | !is.null(TxDb) )

        if( is.null(genomicRegions) ){
            stopifnot( length(bins) == 3 )
            genomicRegions = extractGenomicFeatures(TxDb, upstream_width=bins[1], body_width=bins[2], downstream_width=bins[3])
        }

        GenomeInfoDb::seqlevelsStyle(sampleObject) = "UCSC"
        sampleObject = sortSeqlevels(sampleObject)
        if( any(is.na(seqlengths(sampleObject))) | any(!(seqlevels(GenomeInfoDb::keepStandardChromosomes(TxDb)) %in% seqlevels(sampleObject))) ){
            sampleObject = GenomeInfoDb::keepStandardChromosomes(sampleObject, pruning.mode="coarse")
            if(  length(seqlevels(sampleObject)) < length(seqlevels(GenomeInfoDb::keepStandardChromosomes(TxDb))) )
                seqlevels(sampleObject) = seqlevels(GenomeInfoDb::keepStandardChromosomes(TxDb))
            seqinfo(sampleObject) = seqinfo(GenomeInfoDb::keepStandardChromosomes(TxDb))
        }

        profiles = lapply(seq_along(genomicRegions),
            function(i){
                region_name = names(genomicRegions[i])
                region = genomicRegions[[i]]
                nbin = bins[i]

                tmp = subsetByOverlaps(sampleObject, region, ignore.strand=ignoreStrand)

                weightVec = which(colnames(mcols(tmp))==weightCol)
                if( length(weightVec) > 0 )
                    tmp = rep(tmp, mcols(tmp)[, weightVec])

                tmp_list = lapply(setNames(c("+", "-"), c("Plus", "Minus")),
                    function(std){
                        if( any(strand(tmp)==std) & any(any(strand(region)==std)) ){
                            tmp_std = coverage(tmp[strand(tmp)==std,])

                            # Sort regions (crucial for '-' strand)
                            region_std = sort(unlist(region[strand(region)==std,]))
                            region_std = GRangesList(split(region_std, names(region_std)))
                            if( verbose )
                                message(paste(region_name, "regions on", std, " strand:", length(region_std)))
                            region_std = unlist(region_std, use.names=FALSE)
                            region_std$gene_id = names(region_std)

                            region_dt = data.table(gene_id = region_std$gene_id)
                            region_dt[, exon_index := .I]
                            region_dt[, gene_index := .GRP, by="gene_id"]
                            setkey(region_dt, exon_index)

                            tmp_std = as.data.table(tmp_std[region_std], key="group")
                            if( std=="+" ){
                                tmp_std[, pos := 1:.N, by="group"]
                            }else{
                                tmp_std[, pos := .N:1, by="group"]
                            }
                            setkey(tmp_std, group)

                            tmp_std = tmp_std[region_dt, nomatch=0]
                            tmp_std = tmp_std[order(gene_index, group, pos)]
                            tmp_std[, pos := 1:.N, by="gene_index"]
                            tmp_std[, bin := findInterval(pos, seq(0.5, max(pos)+0.5, length.out=nbin+1)), by="gene_index"]

                            tmp_std[, value := as.numeric(value)]
                            tmp_std[, value := value/(max(pos)/nbin), by="gene_id"] # Normalise by gene bin width
                            tmp_std = tmp_std[, list(value = sum(value)), by=c("gene_id", "bin")]
                           
                           #  tmp_std[, norm := value/sum(value), by="gene_id"] # Normalise by gene total count
                           #  tmp_std = tmp_std[!is.na(norm),]
                        }
                        return(tmp_std)
                    })

                tmp_list = rbindlist(tmp_list, idcol="strand")

                tmp_list[, bin := as.numeric(bin) + sum(bins[0:(i-1)])]

                return(tmp_list)
            })

        names(profiles) = names(genomicRegions)
        profiles = rbindlist(profiles, idcol="region_id")

        profiles = profiles[, list(region_id, strand, bin, value = value/sum(value)), by="gene_id"]

        if( !dropEmpty )
            profiles[!is.finite(value), value := 0]
        profiles = profiles[is.finite(value),]

        profiles[, region_id := factor(region_id, levels=names(genomicRegions))]

        if( collapse )
          profiles = profiles[, list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin")]

        return(profiles)
}
