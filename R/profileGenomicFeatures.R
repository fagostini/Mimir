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
#' @param normType A character string indicating which region normalising method to use. One of 'density' (default), 'max', 'none': can be abbreviated.
#'      Depending on the chosen method the values of each region are normalised using the sum ('density'), the maximum ('max') of the values across the region, or not normalised at all ('none').
#' @param collapse When set to 'TRUE', the profiles are collapsed into a single profile.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
#' @return A \code{\link[data.table:data.table-class]{data.table}} of the normalised binned coverage across the genomic features. Column names are determined by \code{collapse}.
#'
#' @import IRanges
#' @import BiocGenerics
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlevelsStyle keepStandardChromosomes 
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
#' ggplot(profile, aes(x=bin, y=Mean, colour=region_id)) + 
#'   geom_line() +
#'   geom_vline(xintercept=c(100.5, 200.5), linetype="dashed", colour="grey30") +
#'   scale_x_continuous("Relative position",
#'        breaks=c(1, 100.5, 200.5, 300), label=c("-1000", "TSS", "TES", "1000")) +
#'   scale_y_continuous("Average normalised signal") +
#'   coord_cartesian(xlim=c(0, 300)) +
#'   theme_bw() +
#'   theme(legend.position=c(0.9, 0.8), legend.background=element_blank()) +
#'   guides(colour=guide_legend(title=""))

profileGenomicFeatures <- function(genomicRegions=NULL, sampleObject=NULL, TxDb=NULL,
    bins=c(100, 100, 100, 100), weightCol=NULL, ignoreStrand=FALSE, dropEmpty=TRUE,
    normType=c("density", "max", "none"), collapse=TRUE, verbose=TRUE){

        stopifnot( !is.null(genomicRegions) | !is.null(TxDb) )
        stopifnot( !is.null(genomicRegions) & length(genomicRegions)<=length(bins) ) # Early call
        stopifnot( !is.null(weightCol)  & weightCol%in%colnames(mcols(sampleObject)) )
        stopifnot( !any(is.na(seqlengths(sampleObject))) | !is.null(TxDb) )

        normType = match.arg(normType)

        if( is.null(genomicRegions) ){
            stopifnot( length(bins) == 3 )
            genomicRegions = extractGenomicFeatures(TxDb, upstream_width=bins[1], body_width=bins[2], downstream_width=bins[3])
        }

        seqlevelsStyle(sampleObject) = "UCSC"
        sampleObject = sortSeqlevels(sampleObject)
        if( any(is.na(seqlengths(sampleObject))) | any(!(seqlevels(keepStandardChromosomes(TxDb)) %in% seqlevels(sampleObject))) ){
            sampleObject = keepStandardChromosomes(sampleObject, pruning.mode="coarse")
            if(  length(seqlevels(sampleObject)) < length(seqlevels(keepStandardChromosomes(TxDb))) )
                seqlevels(sampleObject) = seqlevels(keepStandardChromosomes(TxDb))
            seqinfo(sampleObject) = seqinfo(keepStandardChromosomes(TxDb))
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

                            if( ignoreStrand ){
                                tmp_std = coverage(tmp)
                            }else{
                                tmp_std = coverage(tmp[strand(tmp)==std,])
                            }
                            # Sort regions (crucial for '-' strand)
                            region_std = sort(unlist(region[strand(region)==std,]))
                            region_std = GRangesList(split(region_std, names(region_std)))
                            if( verbose )
                                message(paste(region_name, "regions on", std, "strand:", length(region_std)))

                            tmp_reg = lapply(seq(1, length(region_std), 1e3),
                                function(rg) {
                                    rg = as.numeric(as.character(rg))
                                    sub = region_std[rg:min(rg + 1e3 - 1, length(region_std)),]
                            
                                    tmp_sub = profileRegions(sub, tmp_std, nbin=nbin)
                            
                                    setnames(tmp_sub, "region_id", "gene_id")

                                    return(tmp_sub)
                                })
                        
                            tmp_reg = rbindlist(tmp_reg)

                            return(tmp_reg)
                        }else{
                            return(NULL)
                        }
                    })

                tmp_list = rbindlist(tmp_list, idcol="strand")

                tmp_list[, bin := as.numeric(bin) + sum(bins[0:(i-1)])]

                return(tmp_list)
            })

        names(profiles) = names(genomicRegions)
        profiles = rbindlist(profiles, idcol="region_id")

        profiles = switch(normType,
            density = profiles[, list(region_id, strand, bin, value = value/sum(value)), by="gene_id"],
            max = profiles[, list(region_id, strand, bin, value = value/max(value)), by="gene_id"],
            none = profiles[, list(region_id, strand, bin, value = value), by="gene_id"])

        if( !dropEmpty )
            profiles[!is.finite(value), value := 0]
        profiles = profiles[is.finite(value),]

        profiles[, region_id := factor(region_id, levels=names(genomicRegions))]

        if( collapse )
          profiles = profiles[, list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin")]

        return(profiles)
}
