#' Calculate metadata profiles across fixed-width features
#' 
#' This function calculates metadata profiles across fixed-width features using a procedure similar to that used in (https://www.sciencedirect.com/science/article/pii/S1097276519303533?via%3Dihub#fig1)[Viphakone et al., 2019].
#' @param fixedRegions A named \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} or list of \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} objects. The \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} width must be the same within the same group, but it can vary between groups.
#' @param sampleObject A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} or \code{\link[GenomicAlignments:GAlignments-class]{GenomicAlignments}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. Required if \code{fixedRegions} or \code{sampleObject} have missing \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param bins An ordered integer vector (must be greater equal that the length of \code{fixedRegions}). The vector order determines the number of bins for each region. If more bins than regions are provided, the additional will be ignored.
#' @param weightCol A single character string. This must be the name of an integer column in the \code{sampleObject} object.
#' @param ignoreStrand When set to 'TRUE', the strand information in \code{sampleObject} is ignored. This does not affect the features in \code{fixedRegions}.
#' @param dropEmpty When set to 'TRUE', the transcripts with no signal in any of their sub-regions will be discarded. When set to 'FALSE', all values of these regions will be set to 0.
#' @param normType A character string indicating which region normalising method to use. One of 'density' (default), 'max', 'none': can be abbreviated.
#'      Depending on the chosen method the values of each region are normalised using the sum ('density'), the maximum ('max') of the values across the region, or not normalised at all ('none').
#' @param collapse When set to 'TRUE', the profiles are collapsed into a single profile.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
#' @return A \code{\link[data.table:data.table-class]{data.table}} of the normalised binned coverage across the fixed features. Column names are determined by \code{collapse}.
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
#' query = genes(TxDb, filter=list(tx_chrom = "chr4"))
#' query = GRangesList(Promoter = promoters(query, upstream=500, downstream=1000))
#' 
#' library("pasillaBamSubset")
#' library("GenomicAlignments")
#' 
#' fl1 <- untreated1_chr4()
#' subject = readGAlignments(fl1)
#'
#' profile = profileFixedFeatures(fixedRegions=query, sampleObject=subject, TxDb=TxDb)
#'
#' library("ggplot2")
#'
#' ggplot(profile, aes(x=bin, y=Mean)) + 
#'   geom_line() +
#'   geom_vline(xintercept=50, linetype="dashed", colour="grey30") +
#'   scale_x_continuous("Relative position",
#'        breaks=c(1, 50.5, 150), label=c("-500", "TSS", "1000")) +
#'   scale_y_continuous("Average normalised signal") +
#'   coord_cartesian(xlim=c(0, 150)) +
#'   theme_bw()

profileFixedFeatures <- function(fixedRegions=NULL, sampleObject=NULL, TxDb=NULL,
    bins=c(150), weightCol=NULL, ignoreStrand=FALSE, dropEmpty=TRUE,
    normType=c("density", "max", "none"), collapse=TRUE, verbose=TRUE){

        stopifnot( !is.null(fixedRegions) )
        stopifnot( !is.null(weightCol)  & weightCol%in%colnames(mcols(sampleObject)) )
        stopifnot( (!any(is.na(seqlengths(fixedRegions))) & !any(is.na(seqlengths(sampleObject)))) & !is.null(TxDb) )
        if( class(fixedRegions)[1]%in%c("CompressedGRangesList", "GRangesList") ){
           stopifnot( length(fixedRegions) <= bins )
           for( i in seq_along(fixedRegions) )
              stopifnot( length(unique(width(fixedRegions[[i]])))==1 )
         }else{
              stopifnot( length(unique(width(fixedRegions)))==1 )
         }

        normType = match.arg(normType)

        seqlevelsStyle(sampleObject) = "UCSC"
        sampleObject = sortSeqlevels(sampleObject)
        if( any(is.na(seqlengths(sampleObject))) ){
            sampleObject = keepStandardChromosomes(sampleObject, pruning.mode="coarse")
            if( length(seqlevels(sampleObject)) < length(seqlevels(keepStandardChromosomes(TxDb))) )
                seqlevels(sampleObject) = seqlevels(keepStandardChromosomes(TxDb))
            seqinfo(sampleObject) = seqinfo(keepStandardChromosomes(TxDb))
        }

        profiles = lapply(seq_along(fixedRegions),
            function(i){
                region_name = names(fixedRegions[i])
                region = fixedRegions[[i]]
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
                            region_std = region[strand(region)==std,]
                            region_std = sort(region_std[sum(width(region_std))>0,])
                            if( verbose )
                                message(paste(region_name, "regions on", std, "strand:", length(region_std)))

                            tmp_reg = lapply(seq(1, length(region_std), 1e3),
                                function(rg) {
                                    rg = as.numeric(as.character(rg))
                                    sub = region_std[rg:min(rg + 1e3 - 1, length(region_std)),]
                            
                                    tmp_sub = profileStrandedRegions(sub, tmp_std, nbin=nbin)
                            
                                    setnames(tmp_sub, "region_id", "feature_id")

                                    return(tmp_sub)
                                })
                        
                            tmp_reg = rbindlist(tmp_reg)

                            return(tmp_reg)
                        }else{
                            return(NULL)
                        }
                    })

                tmp_list = rbindlist(tmp_list, idcol="strand")

                tmp_list[, bin := as.numeric(bin)]

                return(tmp_list)
            })

        names(profiles) = names(fixedRegions)
        profiles = rbindlist(profiles, idcol="region_id")

        profiles = switch(normType,
            density = profiles[, list(region_id, strand, bin, value = value/sum(value)), by="feature_id"],
            max = profiles[, list(region_id, strand, bin, value = value/max(value)), by="feature_id"],
            none = profiles[, list(region_id, strand, bin, value = value), by="feature_id"])

        if( !dropEmpty )
            profiles[!is.finite(value), value := 0]
        profiles = profiles[is.finite(value),]

        profiles[, region_id := factor(region_id, levels=names(fixedRegions))]

        if( collapse )
          profiles = profiles[, list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin")]

        return(profiles)
}
