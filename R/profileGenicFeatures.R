#' Calculate metadata profiles across genic features
#' 
#' This function calculates metadata profiles across genic features using a procedure similar to that used in (https://www.sciencedirect.com/science/article/pii/S1097276519303533?via%3Dihub#fig1)[Viphakone et al., 2019].
#' @param genicRegions A named list of \code{\link[GenomicRanges:GRangesList-class]{GenomicRanges}} objects.
#' @param sampleObject A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} or \code{\link[GenomicAlignments:GAlignments-class]{GenomicAlignments}} object.
#' @param TxDb A \code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}} object. Required if \code{genicRegions} or \code{tx2gene} are not provided. It must contain \code{\link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb}} information.
#' @param tx2gene A \code{\linkS4class{data.frame}} object. The first column must be of transcript identifiers, while the second must be of gene identifiers. Additional columns will be discarded.
#' @param bins An ordered integer vector (must be greater equal that the length of \code{genicRegions}). The vector order determines the number of bins for each region. If more bins than regions are provided, the additional will be ignored.  Default order: 5'UTR, CDS, 3'UTR, Intron.
#' @param weightCol A single character string. This must be the name of an integer column in the \code{sampleObject} object.
#' @param ignoreStrand When set to 'TRUE', the strand information in \code{sampleObject} is ignored. This does not affect the features in \code{genicRegions}.
#' @param dropEmpty When set to 'TRUE', the transcripts with no signal in any of their sub-regions will be discarded. When set to 'FALSE', all values of these regions will be set to 0.
#' @param normType A character string indicating which region normalising method to use. One of 'density' (default), 'max', 'none': can be abbreviated.
#'      Depending on the chosen method the values of each region are normalised using the sum ('density'), the maximum ('max') of the values across the region, or not normalised at all ('none').
#' @param collapseBy A character string indicating which summarising method to use. One of "region" (default), "gene", "transcript", "recursive": can be abbreviated.
#'      Profiles are always calculated per transcript (tx_id), and reported as such when this parameter is set to 'transcript'.
#'      When set to 'gene', profiles are collapsed (using sum, mean and sd) by gene_id, while with 'region' they are collapsed in the same way but by region_id.
#'      When set to 'recursive', profiles are first collapse by gene_id and then by region_id, and sum, mean (pooled) and sd (pooled) are reported. 
#'      Please, be aware that each method will result in a different output, as the number and names of the columns by the chosen method.
#' @param verbose When set to 'TRUE', the function prints diagnostic messages.
#' @return A \code{\link[data.table:data.table-class]{data.table}} of the normalised binned coverage across the genic features. Names and number of columns are determined by \code{collapseBy}.
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
#' query = extractGenicFeatures(TxDb)
#' 
#' library("pasillaBamSubset")
#' library("GenomicAlignments")
#' 
#' fl1 <- untreated1_chr4()
#' subject = readGAlignments(fl1)
#'
#' profile = profileGenicFeatures(genicRegions=query, sampleObject=subject, TxDb=TxDb)
#'
#' library("ggplot2")
#'
#' ggplot(profile, aes(x=bin, y=Mean, colour=region_id)) + 
#'   geom_line() +
#'   geom_vline(xintercept=c(20.5, 120.5), linetype="dashed", colour="grey30") +
#'   scale_x_continuous("Relative position",
#'        breaks=c(10.5, 70.5, 155.5), label=c("5'-UTR", "CDS", "3'-UTR")) +
#'   scale_y_continuous("Average normalised signal") +
#'   coord_cartesian(xlim=c(0, 190)) +
#'   theme_bw() +
#'   theme(legend.position=c(0.9, 0.8), legend.background=element_blank()) +
#'   guides(colour=guide_legend(title=""))

profileGenicFeatures <- function(genicRegions=NULL, sampleObject=NULL, TxDb=NULL, tx2gene=NULL,
    bins=c(20, 100, 70, 100), weightCol=NULL, ignoreStrand=FALSE, dropEmpty=TRUE, normType=c("density", "max", "none"),
    collapseBy=c("region", "gene", "transcript", "recursive"), verbose=TRUE){

        stopifnot( !is.null(genicRegions) | !is.null(TxDb) )
        stopifnot( !is.null(genicRegions) & length(genicRegions)<=length(bins) ) # Early call
        stopifnot( !is.null(tx2gene) | !is.null(TxDb) )
        stopifnot( !is.null(weightCol)  & weightCol%in%colnames(mcols(sampleObject)) )
        stopifnot( !any(is.na(seqlengths(sampleObject))) | !is.null(TxDb) )

        normType = match.arg(normType)
        collapseBy = match.arg(collapseBy)

        if( is.null(tx2gene) ){
            tx2gene = unlist(transcriptsBy(TxDb, by="gene"), use.names=TRUE)
            tx2gene = data.table(tx_id = tx2gene$tx_name, gene_id = names(tx2gene))
        }
        tx2gene = data.table(tx2gene[, 1:2])
        setnames(tx2gene, c("tx_id", "gene_id"))

        if( is.null(genicRegions) ){
            genicRegions = extractGenicFeatures(TxDb, tx2gene=tx2gene, bins=bins)
        }
        stopifnot( !is.null(genicRegions) & length(genicRegions)<=length(bins) ) # Late call

        seqlevelsStyle(sampleObject) = "UCSC"
        sampleObject = sortSeqlevels(sampleObject)
        if( any(is.na(seqlengths(sampleObject))) | any(!(seqlevels(keepStandardChromosomes(TxDb)) %in% seqlevels(sampleObject))) ){
            sampleObject = keepStandardChromosomes(sampleObject, pruning.mode="coarse")
            if(  length(seqlevels(sampleObject)) < length(seqlevels(keepStandardChromosomes(TxDb))) )
                seqlevels(sampleObject) = seqlevels(keepStandardChromosomes(TxDb))
            seqinfo(sampleObject) = seqinfo(keepStandardChromosomes(TxDb))
        }

        profiles = lapply(seq_along(genicRegions),
            function(i){
                region_name = names(genicRegions[i])
                region = genicRegions[[i]]
                nbin = bins[i]

                tmp = subsetByOverlaps(sampleObject, region, ignore.strand=ignoreStrand)

                weightVec = which(colnames(mcols(tmp))==weightCol)
                if( length(weightVec) > 0 )
                    tmp = rep(tmp, mcols(tmp)[, weightVec])

                tmp_list = lapply(setNames(c("+", "-"), c("Plus", "Minus")),
                    function(std){
                        if( any(any(strand(region)==std)) ){
                            
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
                                    
                                    setnames(tmp_sub, "region_id", "tx_id")
                                    setkey(tmp_sub, tx_id)
                                    setkey(tx2gene, tx_id)
                                    tmp_sub[tx2gene, gene_id := gene_id]

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

        names(profiles) = names(genicRegions)
        profiles = rbindlist(profiles, idcol="region_id")

        profiles = switch(normType,
            density = profiles[, list(region_id, strand, bin, value = value/sum(value)), by=c("gene_id", "tx_id")],
            max = profiles[, list(region_id, strand, bin, value = value/max(value)), by=c("gene_id", "tx_id")],
            none = profiles[, list(region_id, strand, bin, value = value), by=c("gene_id", "tx_id")])

        if( !dropEmpty )
            profiles[!is.finite(value), value := 0]
        profiles = profiles[is.finite(value),]

        profiles[, region_id := factor(region_id, levels=names(genicRegions))]

        switch(collapseBy,
             region = profiles[,
                            list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin")],
            gene = profiles[,
                            list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin", "gene_id")],
            transcript = profiles,
            recursive = (profiles[,
                            list(Sum = sum(value), Mean = mean(value), Sd = sd(value)), by=c("region_id", "bin", "gene_id")][,
                            list(Sum = sum(Mean),   PMean = mean(Mean), PSd = sqrt(sum((Mean-mean(Mean)^2))/(.N-1))), by=c("region_id", "bin")]))
}
